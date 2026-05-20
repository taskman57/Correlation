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
        adcI_i      : in STD_LOGIC_VECTOR (15 downto 0);    -- cartesian X
        adcq_i      : in STD_LOGIC_VECTOR (15 downto 0);    -- cartesian Y
        dout_val_o  : out std_logic;
        mag_o       : out STD_LOGIC_VECTOR (15 downto 0);
        phs_o       : out STD_LOGIC_VECTOR (15 downto 0)
    );
end vectoring;

architecture Behavioral of vectoring is

    COMPONENT IQ2Vector
        PORT (
            aclk                        : IN STD_LOGIC;
            aclken                      : IN STD_LOGIC;
            aresetn                     : IN STD_LOGIC;
            s_axis_cartesian_tvalid     : IN STD_LOGIC;
            s_axis_cartesian_tdata      : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
            m_axis_dout_tvalid          : OUT STD_LOGIC;
            m_axis_dout_tdata           : OUT STD_LOGIC_VECTOR(31 DOWNTO 0) 
        );
    END COMPONENT;
    signal dinp_dat_s   : std_logic_vector(31 downto 0);
    signal dout_dat_s   : std_logic_vector(31 downto 0);

begin

IQ2Vector_INST : IQ2Vector
    PORT MAP (
        aclk                        => sys_clk_i,
        aclken                      => '1',
        aresetn                     => not sys_rst_i,
        s_axis_cartesian_tvalid     => inp_val_i,
        s_axis_cartesian_tdata      => dinp_dat_s,
        m_axis_dout_tvalid          => dout_val_o,
        m_axis_dout_tdata           => dout_dat_s
    );
    dinp_dat_s  <= adcq_i & adcI_i;
    mag_o       <= dout_dat_s(15 downto 00);
    phs_o       <= dout_dat_s(31 downto 16);

end Behavioral;
