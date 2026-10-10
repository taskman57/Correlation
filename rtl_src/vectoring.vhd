library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;
library work;
use work.radar_pkg.all;

entity vectoring is
    Port (
        sys_clk_i   : in STD_LOGIC;
        sys_rst_i   : in STD_LOGIC;
        inp_val_i   : in STD_LOGIC; 
        adcI_i      : in STD_LOGIC_VECTOR (18 downto 0);    -- cartesian X
        adcq_i      : in STD_LOGIC_VECTOR (18 downto 0);    -- cartesian Y
        dout_val_o  : out std_logic;
        ampl_o      : out STD_LOGIC_VECTOR (15 downto 0);
        phase_o     : out STD_LOGIC_VECTOR (15 downto 0)
    );
end vectoring;

architecture Behavioral of vectoring is

    signal dinp_dat_s   : std_logic_vector(47 downto 0);
    signal dout_dat_s   : std_logic_vector(31 downto 0);

begin

IQ2Vector_INST : entity work.IQ2Vector
    PORT MAP (
        aclk                        => sys_clk_i,
        aclken                      => '1',
        aresetn                     => not sys_rst_i,
        s_axis_cartesian_tvalid     => inp_val_i,
        s_axis_cartesian_tdata      => dinp_dat_s,
        m_axis_dout_tvalid          => dout_val_o,
        m_axis_dout_tdata           => dout_dat_s
    );
    dinp_dat_s  <= std_logic_vector(resize(signed(adcq_i),24)) & std_logic_vector(resize(signed(adcI_i),24));
    ampl_o      <= dout_dat_s(15 downto 00);
    phase_o     <= dout_dat_s(31 downto 16);

end Behavioral;
