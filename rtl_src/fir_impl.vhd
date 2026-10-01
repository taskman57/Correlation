library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
Library UNISIM;
use UNISIM.vcomponents.all;

use work.radar_pkg.all;

entity fir_impl is
    Port ( 
        clk_i           : in std_logic;
        rst_i           : in std_logic;
        clk_ena_i       : in std_logic;
        cyc_ctr_i       : in integer range 0 to DSP_FOLD_STAGES_C-1;
        data_i          : in std_logic_vector(ADC_BIT_RES_C-1 downto 0);
        coef_i          : in chrp_rom_t;
        rdclk_i         : in std_logic;
        fir_res_o       : out std_logic_vector(15 downto 0);
        fir_vld_o       : out std_logic
    );
end fir_impl;

architecture Behavioral of fir_impl is

    constant REG_OUT_STAGE_C    : integer := 2;
    attribute shreg_extract : string;
    attribute srl_style     : string;

    signal sfr_reg_s        : sfr_fir_t     := (others => (others => '0'));
    attribute shreg_extract of sfr_reg_s    : signal is "yes";
    attribute srl_style     of sfr_reg_s    : signal is "srl_reg";  -- Force SRL implementation with an output FF (srl_reg) for timing closure

    ---- DSP singls
    signal dsp_mode_s       : dsp_mode_t    := (others => (others => '0'));
    signal alu_mod_s        : std_logic_vector(3 downto 0):= (others => '0');
    signal A_D_mod_s        : std_logic_vector(4 downto 0):= (others => '0');

    signal dsp_ainp_s       : ainp_t    := (others => (others => '0'));
    signal dsp_binp_s       : binp_t    := (others => (others => '0'));
    signal dsp_pcin_s       : pcin_t    := (others => (others => '0'));
    signal dsp_pout_s       : pout_t    := (others => (others => '0'));
    signal prod_res_s       : pcin_t    := (others => (others => '0'));
    signal dsp_dinp_s       : dinp_t    := (others => (others => '0'));

    signal reg_out_s        : std_logic_vector(31 downto 0);
    signal data_rdclk_s     : std_logic_vector(31 downto 0);
    signal fir_vld_s        : std_logic;

begin

    OUT_VLD_PROC: process(clk_i)
        variable dly_ctr_v      : integer := 0;
    begin
        if rising_edge(clk_i) then
            if rst_i = '1' then
                dly_ctr_v           := 0;
                fir_vld_s           <= '0';
            else
                if dly_ctr_v < (FIR_LEN_C/2 + DSP_LATANCY_C + 5*REG_OUT_STAGE_C) - 1 then
                    dly_ctr_v   := dly_ctr_v + 1;
                else
                    fir_vld_s   <= '1';
                end if;
            end if;
        end if;
    end process OUT_VLD_PROC;

    SFR_TAP_CTRL_PROC: process(clk_i)
        variable dly_ctr_v      : integer := 0;
    begin
        if rising_edge(clk_i) then
            if clk_ena_i = '1' then
                -- input data shifting
                sfr_reg_s       <= sfr_reg_s(FIR_LEN_C-3 downto 0) & signed(data_i);
            end if;
        end if;
    end process SFR_TAP_CTRL_PROC;

    alu_mod_s   <= x"0";                    -- refere to pg.35 of DSP48E1, ALUMODE
    A_D_mod_s   <= '0' & x"5";              -- refere to pg.31, pg.32 of DSP48E1, INMODE

    DSP_FILT_GEN:
    for dsp_stg in 0 to DSP_FOLD_STAGES_C-1 generate
    begin
        DSP0_IF:
        if dsp_stg = 0  generate
            dsp_ainp_s(dsp_stg)     <= std_logic_vector(resize(signed(data_i),DSP_AINP_LEN_C)) when (cyc_ctr_i = 0) else
                                    std_logic_vector(resize(sfr_reg_s(cyc_ctr_i),DSP_AINP_LEN_C));

            dsp_dinp_s(dsp_stg)     <= std_logic_vector(resize(sfr_reg_s((FIR_LEN_C-2)),DSP_DINP_LEN_C)) when (cyc_ctr_i = 0) else
                                    std_logic_vector(resize(sfr_reg_s((FIR_LEN_C-1) - cyc_ctr_i),DSP_DINP_LEN_C));

            dsp_mode_s(dsp_stg)     <= "000" & "01" & "01"  when (cyc_ctr_i = 0 + (DSP_LATANCY_C - 1 - TIMING_COMP_C)) else     -- Z:0, X&Y: M, refere to pg.34 of DSP48E1
                                    "100" & "01" & "01";                                                                        -- Z:P, X&Y: M, refere to pg.34 of DSP48E1

        end generate;

        dsp_binp_s(dsp_stg)     <= std_logic_vector(resize(coef_i((dsp_stg*DSP_FOLD_STAGES_C) + cyc_ctr_i),DSP_BINP_LEN_C));

        DSPn_IF:
        if dsp_stg > 0 generate
            -- Added '+ dsp_stg' to shift data history back by 1 sample per cascade stage
            dsp_ainp_s(dsp_stg)     <= std_logic_vector(resize(sfr_reg_s(((dsp_stg*DSP_FOLD_STAGES_C)-1 + dsp_stg)),DSP_AINP_LEN_C)) when (cyc_ctr_i = 0) else
                                    std_logic_vector(resize(sfr_reg_s(((dsp_stg*DSP_FOLD_STAGES_C)+cyc_ctr_i + dsp_stg)),DSP_AINP_LEN_C));

            dsp_pcin_s(dsp_stg)     <= prod_res_s(dsp_stg-1);

            dsp_mode_s(dsp_stg)     <= "001" & "01" & "01"  when (cyc_ctr_i = (DSP_LATANCY_C - 1 - TIMING_COMP_C)) else     -- Z:PCIN, X&Y: M
                                    "100" & "01" & "01";                                                                    -- Z:P, X&Y: M
                                    
            -- Added '+ dsp_stg' to shift the symmetric pre-adder data history
            dsp_dinp_s(dsp_stg)     <= std_logic_vector(resize(sfr_reg_s((FIR_LEN_C-2) - ((dsp_stg*DSP_FOLD_STAGES_C)) + dsp_stg),DSP_DINP_LEN_C)) when (cyc_ctr_i = 0) else
                                    std_logic_vector(resize(sfr_reg_s((FIR_LEN_C-1) - ((dsp_stg*DSP_FOLD_STAGES_C) + cyc_ctr_i) + dsp_stg),DSP_DINP_LEN_C));
        end generate;        

        DSP48E1_inst : entity work.DSP_wrapper
        Port map(
            clk_i       => clk_i,
            clk_ena_i   => '1',
            rst_i       => rst_i,
            dsp_mod_i   => dsp_mode_s(dsp_stg),
            alu_mod_i   => alu_mod_s,
            A_D_mod_i   => A_D_mod_s,

            ainp_i      => dsp_ainp_s(dsp_stg),
            binp_i      => dsp_binp_s(dsp_stg),
            pcin_i      => dsp_pcin_s(dsp_stg),
            dinp_i      => dsp_dinp_s(dsp_stg),
            pcout_o     => prod_res_s(dsp_stg),
            pout_o      => dsp_pout_s(dsp_stg)
        );
    end generate;

    reg_fir_proc: process(clk_i)
    begin
        if rising_edge(clk_i) then
            if rst_i = '1' then
                reg_out_s       <= (others => '0');
            else
                if cyc_ctr_i = 3 then
                    reg_out_s   <= dsp_pout_s(DSP_FOLD_STAGES_C-1)(31 downto 0);    -- MSB determined from octave script
                end if;
            end if;
        end if;
    end process reg_fir_proc;

    round_out_proc: process(rdclk_i)
        variable round_fir_v    : std_logic_vector(32 downto 0);
    begin
        if rising_edge(rdclk_i) then
            data_rdclk_s    <= reg_out_s;
            round_fir_v     := conv_round(data_rdclk_s, 16);
            fir_res_o       <= round_fir_v(31 downto 16);
            fir_vld_o       <= fir_vld_s;
        end if;
    end process round_out_proc;

end Behavioral;
