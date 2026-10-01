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

    signal amp_fif_emp_s    : std_logic;
    signal pha_fif_emp_s    : std_logic;
    signal amp_fifo_en_s    : std_logic;
    signal amp_fifo_dat_s   : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal pha_fifo_dat_s   : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal adc_pha_s        : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal fir_amp_s        : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal fir_pha_s        : std_logic_vector(ADC_BIT_RES_C-1 downto 0);

    signal low_lev_s        : std_logic;
    signal mid_lev_s        : std_logic;
    signal hig_lev_s        : std_logic;

    signal cyc_ctr_s        : integer range 0 to DSP_FOLD_STAGES_C - 1:=4;

    signal fir_amp_vld_s    : std_logic;
    signal fir_pha_vld_s    : std_logic;

    -- synthesis translate_off
    file fir_file_amp       : text open write_mode is "../../../../../rtl_sim/matched_filter_amp.dat";
    file fir_file_pha       : text open write_mode is "../../../../../rtl_sim/matched_filter_pha.dat";
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

    adc_amp_inst : entity work.adc_fifo
    PORT MAP (
        rst         => fifo_rst_s,
        wr_clk      => adc_clk_i,
        rd_clk      => dsp_clk_s,
        din         => adc_amp_i,
        wr_en       => adc_vld_i,
        rd_en       => amp_fifo_en_s,
        dout        => amp_fifo_dat_s,
        full        => open,
        empty       => amp_fif_emp_s
    );
    adc_pha_inst : entity work.adc_fifo
    PORT MAP (
        rst         => fifo_rst_s,
        wr_clk      => adc_clk_i,
        rd_clk      => dsp_clk_s,
        din         => adc_pha_i,
        wr_en       => adc_vld_i,
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

    fir_amp_inst: entity work.fir_impl
    Port map( 
        clk_i           => dsp_clk_s,
        rdclk_i         => sys_clk_s,
        rst_i           => fir_rst_s,
        clk_ena_i       => amp_fifo_en_s,
        cyc_ctr_i       => cyc_ctr_s,
        data_i          => amp_fifo_dat_s,
        coef_i          => chrp_ampl_c,
        fir_res_o       => fir_amp_s,
        fir_vld_o       => fir_amp_vld_s
    );

    fir_pha_inst: entity work.fir_impl
    Port map( 
        clk_i           => dsp_clk_s,
        rdclk_i         => sys_clk_s,
        rst_i           => fir_rst_s,
        clk_ena_i       => amp_fifo_en_s,
        cyc_ctr_i       => cyc_ctr_s,
        data_i          => pha_fifo_dat_s,
        coef_i          => chrp_phs_c,
        fir_res_o       => fir_pha_s,
        fir_vld_o       => fir_pha_vld_s
    );

    obj_det_o   <= fir_amp_vld_s and fir_pha_vld_s;

-- synthesis translate_off
    stimulus: process(sys_clk_s)
        variable lin_amp_v  : line;
        variable lin_pha_v  : line;
    begin
        if rising_edge(sys_clk_s) then
            if fir_amp_vld_s = '1' then
                report "Writing amplitude FIR result";
                hwrite(lin_amp_v, fir_amp_s);
                writeline(fir_file_amp, lin_amp_v);
            end if;
            if fir_pha_vld_s = '1' then
                report "Writing phase FIR result";
                hwrite(lin_pha_v, fir_pha_s);
                writeline(fir_file_pha, lin_pha_v);
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
                if unsigned(fir_amp_s) < 2**13-1 then
                    low_lev_s   <= '1';
                    mid_lev_s   <= '0';
                    hig_lev_s   <= '0';
                elsif unsigned(fir_amp_s) > 2**13-1 and unsigned(fir_pha_s) < 2**14-1 then
                    low_lev_s   <= '0';
                    mid_lev_s   <= '1';
                    hig_lev_s   <= '0';
                elsif unsigned(fir_pha_s) > 2**14-1 then
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
