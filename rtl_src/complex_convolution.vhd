library IEEE;
use IEEE.std_logic_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
library work;
use work.radar_pkg.all;
use STD.textio.all;
use ieee.std_logic_textio.all;

entity complex_convolution is
    Port (
        clk_i           : in std_logic;
        rst_i           : in std_logic;
        rd_clk_i        : in std_logic;
        rd_rst_i        : in std_logic;
        clk_ena_i       : in std_logic;
        cyc_ctr_i       : in integer range 0 to DSP_FOLD_STAGES_C-1;
        adc_real_i      : in std_logic_vector (15 downto 0);
        adc_imag_i      : in std_logic_vector (15 downto 0);

        conv_real_o     : buffer std_logic_vector (15 downto 0);
        conv_imag_o     : buffer std_logic_vector (15 downto 0);
        conv_vld_o      : buffer std_logic
    );
end complex_convolution;

architecture rtl of complex_convolution is

    -- defined for complex convolution: (a+jb)*(c+jd)
    signal fir_II_s         : std_logic_vector(31 downto 0);
    signal fir_IQ_s         : std_logic_vector(31 downto 0);
    signal fir_QQ_s         : std_logic_vector(31 downto 0);
    signal fir_QI_s         : std_logic_vector(31 downto 0);

    signal fir_II_ctr_s     : std_logic;
    signal fir_IQ_ctr_s     : std_logic;
    signal fir_QQ_ctr_s     : std_logic;
    signal fir_QI_ctr_s     : std_logic;

    signal cmp_conv_amp_s   : std_logic_vector(32 downto 0):=(others => '0');
    signal cmp_conv_img_s   : std_logic_vector(32 downto 0):=(others => '0');
    signal fir_amp_long_s   : std_logic_vector(32 downto 0):=(others => '0');
    signal fir_pha_long_s   : std_logic_vector(32 downto 0):=(others => '0');

    signal fir_amp_ctr_s    : std_logic;
    signal fir_pha_ctr_s    : std_logic;
    signal conv_ctr_s       : std_logic_vector(3 downto 0);

    signal fir_ampl_s       : std_logic_vector(ADC_BIT_RES_C-1 downto 0);
    signal fir_phas_s       : std_logic_vector(ADC_BIT_RES_C-1 downto 0);

    signal ctr_s            : integer range 0 to 2500-1:=0;

    -- synthesis translate_off
    file cmpx_conv_file     : text open write_mode is "../../../../../rtl_sim/Output_module_quant_complex_convolution.dat";
    -- synthesis translate_on
begin

    -- complex convolution implementation
    fir_II_inst: entity work.fir_impl
    Port map( 
        clk_i           => clk_i,
        rd_clk_i        => rd_clk_i,
        rst_i           => rst_i,
        clk_ena_i       => clk_ena_i,
        cyc_ctr_i       => cyc_ctr_i,
        data_i          => adc_real_i,
        coef_i          => chrp_ampl_c,
        fir_res_o       => fir_II_s,
        fir_vld_o       => fir_II_ctr_s
    );

    fir_IQ_inst: entity work.fir_impl
    Port map( 
        clk_i           => clk_i,
        rd_clk_i        => rd_clk_i,
        rst_i           => rst_i,
        clk_ena_i       => clk_ena_i,
        cyc_ctr_i       => cyc_ctr_i,
        data_i          => adc_real_i,
        coef_i          => chrp_phas_c,
        fir_res_o       => fir_IQ_s,
        fir_vld_o       => fir_IQ_ctr_s
    );

    fir_QQ_inst: entity work.fir_impl
    Port map( 
        clk_i           => clk_i,
        rd_clk_i        => rd_clk_i,
        rst_i           => rst_i,
        clk_ena_i       => clk_ena_i,
        cyc_ctr_i       => cyc_ctr_i,
        data_i          => adc_imag_i,
        coef_i          => chrp_phas_c,
        fir_res_o       => fir_QQ_s,
        fir_vld_o       => fir_QQ_ctr_s
    );

    fir_QI_inst: entity work.fir_impl
    Port map( 
        clk_i           => clk_i,
        rd_clk_i        => rd_clk_i,
        rst_i           => rst_i,
        clk_ena_i       => clk_ena_i,
        cyc_ctr_i       => cyc_ctr_i,
        data_i          => adc_imag_i,
        coef_i          => chrp_ampl_c,
        fir_res_o       => fir_QI_s,
        fir_vld_o       => fir_QI_ctr_s
    );

    complex_calc_proc: process(rd_clk_i)
    begin
        if rising_edge(rd_clk_i) then

            conv_ctr_s     <= conv_ctr_s(conv_ctr_s'left-1 downto 0) & fir_II_ctr_s;
            -- complex convolution real and imaginary part
            cmp_conv_amp_s  <= std_logic_vector(resize(signed(fir_II_s), fir_II_s'length+1) - resize(signed(fir_QQ_s), fir_QQ_s'length+1));
            cmp_conv_img_s  <= std_logic_vector(resize(signed(fir_IQ_s), fir_IQ_s'length+1) + resize(signed(fir_QI_s), fir_QI_s'length+1));

            -- convergent rounding to the nearest even
            fir_amp_long_s  <= conv_round(std_logic_vector(signed(cmp_conv_amp_s(cmp_conv_amp_s'left-1 downto 0))), 16);
            fir_pha_long_s  <= conv_round(std_logic_vector(signed(cmp_conv_img_s(cmp_conv_img_s'left-1 downto 0))), 16);

            fir_ampl_s      <= fir_amp_long_s(31 downto 16);
            fir_phas_s      <= fir_pha_long_s(31 downto 16);
        end if;
    end process;

    process(rd_clk_i)
        variable cnt_en_v   : boolean := false;
    begin
        if rising_edge(rd_clk_i) then
            if rd_rst_i = '1' then
                cnt_en_v    := false;
                ctr_s       <= 0;
                conv_vld_o  <= '0';
                conv_real_o <= (others => '0');
                conv_imag_o <= (others => '0');
            else
                if conv_ctr_s(2) = '1' and conv_ctr_s(3) = '0' then     -- just once after start
                    cnt_en_v   := true;
                end if;
                if cnt_en_v then
                    if ctr_s < 2500-1 then
                        ctr_s   <= ctr_s + 1;
                    else
                        ctr_s   <= 0;
                    end if;
                    if ctr_s > 50-1 then
                        conv_vld_o      <= '1';
                        conv_real_o     <= fir_ampl_s;
                        conv_imag_o     <= fir_phas_s;
                    else
                        conv_vld_o      <= '0';
                    end if;
                end if;
            end if;
        end if;
    end process;

-- synthesis translate_off
    simulation: process(rd_clk_i)
        variable lin_compx_v  : line;
    begin
        if rising_edge(rd_clk_i) then
            -- if conv_ctr_s(2) = '1' then
            if conv_vld_o = '1' then
                report "Writing amplitude FIR result";
                -- write(lin_compx_v, to_hex(cmp_conv_amp_s));
                -- write(lin_compx_v, to_integer(signed(cmp_conv_amp_s)));
                -- write(lin_compx_v, to_integer(signed(fir_ampl_s)));
                write(lin_compx_v, to_integer(signed(conv_real_o)));

                write(lin_compx_v, string'(" "));

                -- write(lin_compx_v, to_hex(cmp_conv_img_s));
                -- write(lin_compx_v, to_integer(signed(cmp_conv_img_s)));
                -- write(lin_compx_v, to_integer(signed(fir_phas_s)));
                write(lin_compx_v, to_integer(signed(conv_imag_o)));

                writeline(cmpx_conv_file, lin_compx_v);
            end if;
        end if;
    end process;
-- synthesis translate_on

end rtl;
