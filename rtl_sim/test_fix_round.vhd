library IEEE;
use IEEE.STD_LOGIC_1164.ALL;
use IEEE.NUMERIC_STD.ALL;

entity test_fix_round is
    generic(
        inp_int_g   : integer :=3;
        inp_frc_g   : integer :=0
    );
end test_fix_round;

architecture Behavioral of test_fix_round is
component fix_round is
    generic(
        inp_int_g   : integer :=3;
        inp_frc_g   : integer :=0
    );
    port(
        clk_i   : in std_logic;
        numb_i  : in std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
        numb_o  : out std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
        ovf_o   : out std_logic
    );
end component;

    ---- input
    signal clk_i   : std_logic:='0';
    signal numb_i  : std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0):=(others => '0');

    ---- output
    signal numb_o  : std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
    signal ovf_o   : std_logic;

    ---- constants
    constant clk_freq   : time:= 10 ns;
    type    sample_t    is array(0 to 7) of std_logic_vector((inp_int_g+inp_frc_g)-1 downto 0);
    -- constant test_fix_c : sample_t :=(
        -- "00100",    -- 1.0
        -- "00101",    -- 1.25
        -- "00110",    -- 1.5
        -- "00111",    -- 1.75
        -- "01000",    -- 2.0
        -- "01010",    -- 2.5
        -- "01101",    -- 3.25
        -- "01010",    -- 2.5
        -- "01111",    -- 3.75
        -- "10010",    -- -3.5
        -- "00011",    -- 0.75
        -- "11101",    -- -0.75
        -- "11010",    -- -1.5
        -- "01110",    -- 3.5
        -- "10111",    -- -2.5
        -- "10000"     -- -4
    -- );
    constant test_fix_c : sample_t :=(
        "001",    -- 1
        "010",    -- 2
        "011",    -- 3
        "000",    -- 0
        "111",    -- -1
        "110",    -- -2
        "101",    -- -3
        "100"     -- -4
    );

begin

    clk_itimulu:
    process
    begin
        clk_i <= not clk_i;
        wait for clk_freq/2;
    end process;

    inp_stimulu:
    process(clk_i)
        variable ctr_v  : integer range 0 to 15:=0;
        variable pwu_v  : integer range 0 to 15:=0;
    begin
        if rising_edge(clk_i) then
            if pwu_v < 7 then   -- power up delay
                pwu_v   := pwu_v + 1;
                numb_i  <= (others => '0');
            else
                numb_i <= test_fix_c(ctr_v);
                if ctr_v < 7 then
                    ctr_v   := ctr_v +1;
                end if;
            end if;
        end if;
    end process;

    dut:fix_round
        generic map(
            inp_int_g => inp_int_g,
            inp_frc_g => inp_frc_g
        )
        port map(
            clk_i   => clk_i,
            numb_i  => numb_i,
            numb_o  => numb_o,
            ovf_o   => ovf_o
        );
    
end Behavioral;
