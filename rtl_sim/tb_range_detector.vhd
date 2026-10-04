library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use work.radar_pkg.all;
use work.ref_adc_pkg.all;



entity tb_range_detector is
    generic(SIM_LP_FILTER   : boolean := false);
end tb_range_detector;

architecture sim of tb_range_detector is

    -- inputs
    signal clk_p            : std_logic:='1';
    signal clk_n            : std_logic:='0';
    signal rst_i            : std_logic:='1';
    signal adc_clk_i        : std_logic:='0';
    signal adc_vld_i        : std_logic:='0';
    signal adc_amp_i        : std_logic_vector(ADC_BIT_RES_C-1 downto 0):=(others => '0');
    signal adc_pha_i        : std_logic_vector(ADC_BIT_RES_C-1 downto 0):=(others => '0');
    signal pulse_i          : std_logic:='0';

    --  outputs
    signal low_lev_o        : std_logic;
    signal mid_lev_o        : std_logic;
    signal hig_lev_o        : std_logic;
    signal obj_det_o        : std_logic;

    -- constants
    constant clk_per_c      : time := 20 ns;    -- 50MHz system oscilator clock frequency

begin

    dut: entity work.range_detector
    port map(
        clk_p       => clk_p,
        clk_n       => clk_n,
        rst_i       => rst_i,
        adc_clk_i   => adc_clk_i,
        adc_vld_i   => adc_vld_i,
        adc_amp_i   => adc_amp_i,
        adc_pha_i   => adc_pha_i,
        pulse_i     => pulse_i,
        low_lev_o   => low_lev_o,
        mid_lev_o   => mid_lev_o,
        hig_lev_o   => hig_lev_o,
        obj_det_o   => obj_det_o
    );

    clk_stim: process
    begin
        clk_p <= not clk_p;
        wait for clk_per_c/2;
        clk_n <= clk_p;
    end process;

    rst_stim:
        rst_i <= '0' after 5.1*clk_per_c;

    adc_clk_i <= clk_p;

    adc_data_stim: process(adc_clk_i)
        variable stu_dly_v      : integer range 0 to 100-1:=0;
        variable amp_ctr_v      : integer range 0 to 2500-1:=0;
        variable pha_ctr_v      : integer range 0 to 2500-1:=0;
        variable p_num_v        : integer range 0 to adc_real_c'length-1:=0;    -- Np in Correlation script 
    begin
        if rising_edge(adc_clk_i) then
            if rst_i = '1' then
                amp_ctr_v   := 0;
                pha_ctr_v   := 0;
                stu_dly_v   := 0;
                adc_vld_i   <= '0';
                pulse_i     <= '0';
            else
                if stu_dly_v < 99 then
                    stu_dly_v   := stu_dly_v + 1;
                else
                    adc_vld_i   <= '1';
                    if amp_ctr_v = 0 then
                        pulse_i     <= '1';
                    else
                        pulse_i     <= '0';
                    end if;
                    if SIM_LP_FILTER then
                        adc_amp_i       <= std_logic_vector(noisy_dat_c(amp_ctr_v) - 32768);
                        adc_pha_i       <= std_logic_vector(noisy_dat_c(pha_ctr_v) - 32768);
                        if amp_ctr_v < 100-1 then
                            amp_ctr_v   := amp_ctr_v + 1;
                        else
                            amp_ctr_v   := 0;
                        end if;
                        if pha_ctr_v < 100-1 then
                            pha_ctr_v   := pha_ctr_v + 1;
                        else
                            pha_ctr_v   := 0;
                        end if;
                    else
                        adc_amp_i       <= std_logic_vector(adc_real_c(p_num_v)(amp_ctr_v));
                        adc_pha_i       <= std_logic_vector(adc_imag_c(p_num_v)(pha_ctr_v));
                        if amp_ctr_v < 2500-1 then
                            amp_ctr_v       := amp_ctr_v + 1;
                            pha_ctr_v       := pha_ctr_v + 1;
                        else
                            if p_num_v < adc_real_c'length-1 then
                                amp_ctr_v   := 0;
                                pha_ctr_v   := 0;
                                p_num_v     := p_num_v + 1;
                            else
                                adc_vld_i   <= '0';
                            end if;
                        end if;
                    end if;
                end if;
            end if;
        end if;
    end process;

end sim;
