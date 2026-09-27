library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
library UNISIM;
use UNISIM.VComponents.all;
use work.radar_pkg.all;


entity DSP_wrapper is
    Port (
        clk_i       : in std_logic;
        clk_ena_i   : in std_logic;
        rst_i       : in std_logic;
        dsp_mod_i   : in std_logic_vector(6 downto 0);
        alu_mod_i   : in std_logic_vector(3 downto 0);
        A_D_mod_i   : in std_logic_vector(4 downto 0);

        ainp_i      : in std_logic_vector(DSP_AINP_LEN_C-1 downto 0);
        binp_i      : in std_logic_vector(DSP_BINP_LEN_C-1 downto 0);
        pcin_i      : in std_logic_vector(DSP_CINP_LEN_C-1 downto 0);
        dinp_i      : in std_logic_vector(DSP_DINP_LEN_C-1 downto 0);

        pcout_o     : out std_logic_vector(DSP_POUT_LEN_C-1 downto 0);
        pout_o      : out std_logic_vector(DSP_POUT_LEN_C-1 downto 0)
    );
end DSP_wrapper;

architecture RTL of DSP_wrapper is

begin

    DSP48E1_inst : DSP48E1
        generic map (
            -- Feature Control Attributes: Data Path Selection
            A_INPUT => "DIRECT",                  -- Selects A input source, "DIRECT" (A port) or "CASCADE" (ACIN port)
            B_INPUT => "DIRECT",                  -- Selects B input source, "DIRECT" (B port) or "CASCADE" (BCIN port)
            USE_DPORT => TRUE,                   -- Select D port usage (TRUE or FALSE)
            USE_MULT => "MULTIPLY",               -- Select multiplier usage ("MULTIPLY", "DYNAMIC", or "NONE")
            USE_SIMD => "ONE48",                  -- SIMD selection ("ONE48", "TWO24", "FOUR12")
            -- Pattern Detector Attributes: Pattern Detection Configuration
            AUTORESET_PATDET => "NO_RESET",       -- "NO_RESET", "RESET_MATCH", "RESET_NOT_MATCH" 
            MASK => X"3fffffffffff",              -- 48-bit mask value for pattern detect (1=ignore)
            PATTERN => X"000000000000",           -- 48-bit pattern match for pattern detect
            SEL_MASK => "MASK",                   -- "C", "MASK", "ROUNDING_MODE1", "ROUNDING_MODE2" 
            SEL_PATTERN => "PATTERN",             -- Select pattern value ("PATTERN" or "C")
            USE_PATTERN_DETECT => "NO_PATDET",    -- Enable pattern detect ("PATDET" or "NO_PATDET")
            -- Register Control Attributes: Pipeline Register Configuration
            ACASCREG => 1,                        -- Number of pipeline stages between A/ACIN and ACOUT (0, 1 or 2)
            ADREG => 1,                           -- Number of pipeline stages for pre-adder (0 or 1)
            ALUMODEREG => 0,                      -- Number of pipeline stages for ALUMODE (0 or 1)
            AREG => 1,                            -- Number of pipeline stages for A (0, 1 or 2)
            BCASCREG => 1,                        -- Number of pipeline stages between B/BCIN and BCOUT (0, 1 or 2)
            BREG => 2,                            -- Number of pipeline stages for B (0, 1 or 2)
            CARRYINREG => 0,                      -- Number of pipeline stages for CARRYIN (0 or 1)
            CARRYINSELREG => 0,                   -- Number of pipeline stages for CARRYINSEL (0 or 1)
            CREG => 1,                            -- Number of pipeline stages for C (0 or 1)
            DREG => 1,                            -- Number of pipeline stages for D (0 or 1)
            INMODEREG => 0,                       -- Number of pipeline stages for INMODE (0 or 1)
            MREG => 1,                            -- Number of multiplier pipeline stages (0 or 1)
            OPMODEREG => 1,                       -- Number of pipeline stages for OPMODE (0 or 1)
            PREG => 1                             -- Number of pipeline stages for P (0 or 1)
        )
        port map (
            -- Cascade: 30-bit (each) output: Cascade Ports
            ACOUT => open,                        -- 30-bit output: A port cascade output
            BCOUT => open,                        -- 18-bit output: B port cascade output
            CARRYCASCOUT => open,                 -- 1-bit output: Cascade carry output
            MULTSIGNOUT => open,                  -- 1-bit output: Multiplier sign cascade output
            PCOUT => pcout_o,                     -- 48-bit output: Cascade output
            -- Control: 1-bit (each) output: Control Inputs/Status Bits
            OVERFLOW => open,                     -- 1-bit output: Overflow in add/acc output
            PATTERNBDETECT => open,               -- 1-bit output: Pattern bar detect output
            PATTERNDETECT => open,                -- 1-bit output: Pattern detect output
            UNDERFLOW => open,                    -- 1-bit output: Underflow in add/acc output
            -- Data: 4-bit (each) output: Data Ports
            CARRYOUT => open,                     -- 4-bit output: Carry output
            P => pout_o,                          -- 48-bit output: Primary data output
            -- Cascade: 30-bit (each) input: Cascade Ports
            ACIN => (others => '0'),              -- 30-bit input: A cascade data input
            BCIN => (others => '0'),              -- 18-bit input: B cascade input
            CARRYCASCIN => '0',                   -- 1-bit input: Cascade carry input
            MULTSIGNIN => '0',                    -- 1-bit input: Multiplier sign input
            PCIN => pcin_i,                       -- 48-bit input: P cascade input
            -- Control: 4-bit (each) input: Control Inputs/Status Bits
            ALUMODE => alu_mod_i,                 -- 4-bit input: ALU control input
            CARRYINSEL => "000",                  -- 3-bit input: Carry select input
            CLK => clk_i,                         -- 1-bit input: Clock input
            INMODE => A_D_mod_i,                  -- 5-bit input: INMODE control input
            OPMODE => dsp_mod_i,                  -- 7-bit input: Operation mode input
            -- Data: 30-bit (each) input: Data Ports
            A => ainp_i,                          -- 30-bit input: A data input
            B => binp_i,                          -- 18-bit input: B data input
            C => (others => '0'),                 -- 48-bit input: C data input
            CARRYIN => '0',                       -- 1-bit input: Carry input signal
            D => dinp_i,                          -- 25-bit input: D data input
            -- Reset/Clock Enable: 1-bit (each) input: Reset/Clock Enable Inputs
            CEA1 => clk_ena_i,                    -- 1-bit input: Clock enable input for 1st stage AREG
            CEA2 => '0',                          -- 1-bit input: Clock enable input for 2nd stage AREG
            CEAD => clk_ena_i,                    -- 1-bit input: Clock enable input for ADREG
            CEALUMODE => clk_ena_i,               -- 1-bit input: Clock enable input for ALUMODE
            CEB1 => clk_ena_i,                    -- 1-bit input: Clock enable input for 1st stage BREG
            CEB2 => clk_ena_i,                    -- 1-bit input: Clock enable input for 2nd stage BREG
            CEC => clk_ena_i,                     -- 1-bit input: Clock enable input for CREG
            CECARRYIN => '0',                     -- 1-bit input: Clock enable input for CARRYINREG
            CECTRL => clk_ena_i,                  -- 1-bit input: Clock enable input for OPMODEREG and CARRYINSELREG
            CED => clk_ena_i,                     -- 1-bit input: Clock enable input for DREG
            CEINMODE => '1',                      -- 1-bit input: Clock enable input for INMODEREG
            CEM => clk_ena_i,                     -- 1-bit input: Clock enable input for MREG
            CEP => clk_ena_i,                     -- 1-bit input: Clock enable input for PREG
            RSTA => rst_i,                        -- 1-bit input: Reset input for AREG
            RSTALLCARRYIN => rst_i,               -- 1-bit input: Reset input for CARRYINREG
            RSTALUMODE => rst_i,                  -- 1-bit input: Reset input for ALUMODEREG
            RSTB => rst_i,                        -- 1-bit input: Reset input for BREG
            RSTC => rst_i,                        -- 1-bit input: Reset input for CREG
            RSTCTRL => rst_i,                     -- 1-bit input: Reset input for OPMODEREG and CARRYINSELREG
            RSTD => rst_i,                        -- 1-bit input: Reset input for DREG and ADREG
            RSTINMODE => rst_i,                   -- 1-bit input: Reset input for INMODEREG
            RSTM => rst_i,                        -- 1-bit input: Reset input for MREG
            RSTP => rst_i                         -- 1-bit input: Reset input for PREG
        );

end RTL;
