library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
use ieee.math_real.all;

entity IQ2Vect_tb is
--  Port ( );
end IQ2Vect_tb;

architecture Behavioral of IQ2Vect_tb is

component vectoring is
    Port (
        sys_clk_i   : in STD_LOGIC;
        sys_rst_i   : in STD_LOGIC;
        inp_val_i   : in STD_LOGIC; 
        inphase_i   : in STD_LOGIC_VECTOR (15 downto 0);    -- cartesian X
        quadrat_i   : in STD_LOGIC_VECTOR (15 downto 0);    -- cartesian Y
        dout_val_o  : out std_logic;
        mag_o       : out STD_LOGIC_VECTOR (15 downto 0);
        phs_o       : out STD_LOGIC_VECTOR (15 downto 0)
    );
end component;

    -- input
        signal sys_clk_i    : STD_LOGIC := '0';
        signal sys_rst_i    : STD_LOGIC := '1';
        signal inp_val_i    : STD_LOGIC := '0'; 
        signal inphase_i    : STD_LOGIC_VECTOR (15 downto 0) := (others => '0');
        signal quadrat_i    : STD_LOGIC_VECTOR (15 downto 0) := (others => '0');

    -- output
        signal dout_val_o   : std_logic;
        signal mag_o        : STD_LOGIC_VECTOR (15 downto 0);
        signal phs_o        : STD_LOGIC_VECTOR (15 downto 0);

    -- stimulu
        constant clk_per    : time:= 10 ns;

begin

    clk_STIMULU:
    process
    begin
        sys_clk_i <= not sys_clk_i;
        wait for clk_per/2;
    end process;

    sys_rst_i   <= '0' after 1.1*clk_per;

    process(sys_clk_i)
        variable teta_v     : real := 0.0;
        variable cos_v      : real := 0.0;
        variable sin_v      : real := 0.0;
    begin
        if rising_edge(sys_clk_i) then
            if sys_rst_i = '1' then
                inp_val_i   <= '0';
                inphase_i   <= (others => '0');
                quadrat_i   <= (others => '0');
                cos_v       := 0.0;
                sin_v       := 0.0;
                teta_v      := 0.0;
            else
                inp_val_i   <= '1';
                if teta_v <= 1.0 then
                    teta_v := teta_v + 0.001;
                else
                    teta_v := 0.001;
                end if;
                cos_v       := COS(2*MATH_PI*teta_v);
                sin_v       := SIN(2*MATH_PI*teta_v);
                -- Change 15 to 14 so that 1.0 maps to 16384 (the CORDIC's "1.0")
                inphase_i   <= std_logic_vector(to_signed(integer(ceil((2.0**14-0.0) * cos_v)),16));
                quadrat_i   <= std_logic_vector(to_signed(integer(ceil((2.0**14-0.0) * sin_v)),16));
            end if;
        end if;
    end process;

    -- inphase_i   <= "0010000000000000";
    -- quadrat_i   <= "0010000000000000";
    

    DUT: component vectoring
    port map(
        sys_clk_i   => sys_clk_i,
        sys_rst_i   => sys_rst_i,
        inp_val_i   => inp_val_i,
        inphase_i   => inphase_i,
        quadrat_i   => quadrat_i,
        dout_val_o  => dout_val_o,
        mag_o       => mag_o,
        phs_o       => phs_o
    );

end Behavioral;
