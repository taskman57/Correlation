library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

library work;
use work.radar_pkg.all;

Library xpm;
use xpm.vcomponents.all;
use STD.textio.all;
use ieee.std_logic_textio.all;

entity range_detector is
    Port (
        clk_p           : in std_logic;
        clk_n           : in std_logic;
        rst_i           : in std_logic;
        adc_clk_i       : in std_logic;
        adc_vld_i       : in std_logic;
        adc_amp_i       : in std_logic_vector(ADC_BIT_RES_C-1 downto 0);
        adc_pha_i       : in std_logic_vector(ADC_BIT_RES_C-1 downto 0);
        pulse_i         : in std_logic;

        --  processing results
        low_lev_o       : out std_logic;
        mid_lev_o       : out std_logic;
        hig_lev_o       : out std_logic;
        obj_det_o       : out std_logic
    );
end range_detector;

architecture Behavioral of range_detector is

    signal dsp_clk_s        : std_logic;
    signal sys_clk_s        : std_logic;
    signal dsp_rst_s        : std_logic;
    signal fir_rst_s        : std_logic;
    signal sys_rst_s        : std_logic;
    signal clk_lck_s        : std_logic;
    signal fifo_rst_syn_s   : std_logic;
    signal fifo_rst_s       : std_logic;

    signal adc_vld_s        : std_logic;
    signal amp_fif_emp_s    : std_logic;
    signal pha_fif_emp_s    : std_logic;
    signal amp_fifo_en_s    : std_logic;

    signal adc_amp_s        : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal adc_pha_s        : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal amp_fifo_dat_s   : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal pha_fifo_dat_s   : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal conv_real_s      : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal conv_imag_s      : std_logic_vector(ADC_BIT_RES_C-1 downto 0);

    signal acc_conv_real_s  : std_logic_vector(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0);
    signal acc_conv_imag_s  : std_logic_vector(ADC_BIT_RES_C+ACCUM_BIT_GROWTH_C-1 downto 0);

    signal conv_vld_s       : std_logic;
    signal acc_conv_vld_s   : std_logic;
    signal low_lev_s        : std_logic;
    signal mid_lev_s        : std_logic;
    signal hig_lev_s        : std_logic;

    signal cyc_ctr_s        : integer range 0 to DSP_FOLD_STAGES_C - 1:=4;

    -- synthesis translate_off
    file cmpx_conv_file     : text open write_mode is "../../../../../rtl_sim/coherent_summation.dat";
    -- synthesis translate_on

begin

    ----    system clock management
    clk_dsp_inst : entity work.clk_dsp
    port map (
        -- Clock out ports  
        clk_out1        => dsp_clk_s,   -- 250 MHz
        clk_out2        => sys_clk_s,   -- 50 MHz
        -- Status and control signals
        reset           => '0',
        locked          => clk_lck_s,
        -- Clock in ports
        clk_in1_p       => clk_p,
        clk_in1_n       => clk_n
    );

    ----    reset synchronizers
    dsp_async_rst_inst : xpm_cdc_async_rst
    generic map (
        DEST_SYNC_FF    => 9,   -- DECIMAL; range: 2-10
        INIT_SYNC_FF    => 0,   -- DECIMAL; 0=disable simulation init values, 1=enable simulation init values
        RST_ACTIVE_HIGH => 0    -- DECIMAL; 0=active low reset, 1=active high reset
    )
    port map (
        dest_arst       => dsp_rst_s,       -- 1-bit output: src_arst asynchronous reset signal synchronized to destination
                                            -- clock domain. This output is registered. NOTE: Signal asserts asynchronously
                                            -- but deasserts synchronously to dest_clk. Width of the reset signal is at least
                                            -- (DEST_SYNC_FF*dest_clk) period.

        dest_clk        => dsp_clk_s,       -- 1-bit input: Destination clock.
        src_arst        => rst_i            -- 1-bit input: Source asynchronous reset signal.
    );

    sys_async_rst_inst : xpm_cdc_async_rst
    generic map (
        DEST_SYNC_FF    => 9,   -- DECIMAL; range: 2-10
        INIT_SYNC_FF    => 0,   -- DECIMAL; 0=disable simulation init values, 1=enable simulation init values
        RST_ACTIVE_HIGH => 0    -- DECIMAL; 0=active low reset, 1=active high reset
    )
    port map (
        dest_arst       => sys_rst_s,       -- 1-bit output: src_arst asynchronous reset signal synchronized to destination
                                            -- clock domain. This output is registered. NOTE: Signal asserts asynchronously
                                            -- but deasserts synchronously to dest_clk. Width of the reset signal is at least
                                            -- (DEST_SYNC_FF*dest_clk) period.

        dest_clk        => sys_clk_s,       -- 1-bit input: Destination clock.
        src_arst        => rst_i            -- 1-bit input: Source asynchronous reset signal.
    );

    -- RX blanking: switching duplexer between TX/RX, when TX ADC value must be cleared
    rx_blanking_proc: process(adc_clk_i)
        variable blank_ctr_v    : integer range 0 to 50-1 := 0;
    begin
        if rising_edge(adc_clk_i) then
            adc_amp_s   <= adc_amp_i;
            adc_pha_s   <= adc_pha_i;
            if pulse_i = '1' then
                blank_ctr_v := 0;
                adc_amp_s   <= (others => '0');
                adc_pha_s   <= (others => '0');
            elsif blank_ctr_v < 50-1 then
                blank_ctr_v := blank_ctr_v + 1;
                adc_amp_s   <= (others => '0');
                adc_pha_s   <= (others => '0');
            end if;
            adc_vld_s   <= adc_vld_i;
        end if;
    end process;

    adc_real_inst : entity work.adc_fifo
    PORT MAP (
        rst         => fifo_rst_s,
        wr_clk      => adc_clk_i,
        rd_clk      => dsp_clk_s,
        din         => adc_amp_s,
        wr_en       => adc_vld_s,
        rd_en       => amp_fifo_en_s,
        dout        => amp_fifo_dat_s,
        full        => open,
        empty       => amp_fif_emp_s
    );
    adc_imag_inst : entity work.adc_fifo
    PORT MAP (
        rst         => fifo_rst_s,
        wr_clk      => adc_clk_i,
        rd_clk      => dsp_clk_s,
        din         => adc_pha_s,
        wr_en       => adc_vld_s,
        rd_en       => amp_fifo_en_s,
        dout        => pha_fifo_dat_s,
        full        => open,
        empty       => pha_fif_emp_s
    );

    process(sys_clk_s, sys_rst_s)
    begin
        if sys_rst_s = '1' then
            fifo_rst_s      <= '1';
            fifo_rst_syn_s  <= '1';
        elsif rising_edge(sys_clk_s) then
            fifo_rst_syn_s  <= '0';
            fifo_rst_s      <= fifo_rst_syn_s;
        end if;
    end process;

    process(dsp_clk_s)
    begin
        if rising_edge(dsp_clk_s) then
            if dsp_rst_s = '1' then
                amp_fifo_en_s       <= '0';
                cyc_ctr_s           <= 4;
                fir_rst_s           <= '1';
            else
                amp_fifo_en_s       <= '0';
                if amp_fif_emp_s = '0' then
                    fir_rst_s       <= '0';
                    if cyc_ctr_s = 4 then
                        amp_fifo_en_s   <= '1';
                        cyc_ctr_s       <= 0;
                    end if;
                end if;
                if cyc_ctr_s < 4 then
                    cyc_ctr_s       <= cyc_ctr_s + 1;
                end if;
            end if;
        end if;
    end process;

    complex_conv_inst: entity work.complex_convolution
    Port map(
        clk_i           => dsp_clk_s,
        rst_i           => dsp_rst_s,
        rd_clk_i        => sys_clk_s,
        rd_rst_i        => sys_rst_s,
        clk_ena_i       => amp_fifo_en_s, 
        cyc_ctr_i       => cyc_ctr_s, 
        adc_real_i      => amp_fifo_dat_s,
        adc_imag_i      => pha_fifo_dat_s,

        conv_real_o     => conv_real_s,
        conv_imag_o     => conv_imag_s,
        conv_vld_o      => conv_vld_s
    );

    real_coherent_sum_inst: entity work.coherent_sum
    Port map(
        clk_i   => sys_clk_s,
        rst_i   => sys_rst_s,
        data_i  => conv_real_s,
        ena_i   => conv_vld_s,
        sum_o   => acc_conv_real_s,
        vld_o   => acc_conv_vld_s
    );

    imag_coherent_sum_inst: entity work.coherent_sum
    Port map(
        clk_i   => sys_clk_s,
        rst_i   => sys_rst_s,
        data_i  => conv_imag_s,
        ena_i   => conv_vld_s,
        sum_o   => acc_conv_imag_s,
        vld_o   => open
    );
    obj_det_o   <= acc_conv_vld_s;
-- synthesis translate_off
    simulation: process(sys_clk_s)
        variable lin_compx_v  : line;
    begin
        if rising_edge(sys_clk_s) then
            if acc_conv_vld_s = '1' then
                report "Writing coherent summation result";
                write(lin_compx_v, to_integer(signed(acc_conv_real_s)));

                write(lin_compx_v, string'(" "));

                write(lin_compx_v, to_integer(signed(acc_conv_imag_s)));

                writeline(cmpx_conv_file, lin_compx_v);
            end if;
        end if;
    end process;
-- synthesis translate_on

    just_4_test:
    process(sys_clk_s)
    begin
        if rising_edge(sys_clk_s) then
            if sys_rst_s = '1' then
                low_lev_s   <= '0';
                mid_lev_s   <= '0';
                hig_lev_s   <= '0';
            else
                if unsigned(acc_conv_real_s) < 2**13-1 then
                    low_lev_s   <= '1';
                    mid_lev_s   <= '0';
                    hig_lev_s   <= '0';
                elsif unsigned(acc_conv_real_s) > 2**13-1 and unsigned(acc_conv_imag_s) < 2**14-1 then
                    low_lev_s   <= '0';
                    mid_lev_s   <= '1';
                    hig_lev_s   <= '0';
                elsif unsigned(acc_conv_imag_s) > 2**14-1 then
                    low_lev_s   <= '0';
                    mid_lev_s   <= '0';
                    hig_lev_s   <= '1';
                end if;
                low_lev_o       <= low_lev_s;
                mid_lev_o       <= mid_lev_s;
                hig_lev_o       <= hig_lev_s;
            end if;
        end if;
    end process;

end Behavioral;
