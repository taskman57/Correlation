library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
library work;
use work.radar_pkg.all;


entity fix_round is
    generic(
        inp_int_g   : integer :=3;
        inp_frc_g   : integer :=0
    );
    port(
        clk_i   : in std_logic;
        numb_i  : in  std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
        numb_o  : out std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
        ovf_o   : out std_logic
    );
end entity;

architecture rtl of fix_round is

begin

    process(clk_i)
        variable round_v    : std_logic_vector((inp_int_g+inp_frc_g)-0 downto 0):=(others => '0');
    begin
        if rising_edge(clk_i) then
            round_v := conv_round(numb_i, inp_frc_g);
            numb_o  <= round_v((inp_int_g+inp_frc_g)-1 downto 0);
            ovf_o   <= round_v(inp_int_g+inp_frc_g);
        end if;
    end process;

end architecture;