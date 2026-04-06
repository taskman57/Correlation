library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package radar_pkg is

    type chrp_rom_t is array(0 to 49) of signed(15 downto 0);

    -- These Inphase & quadrature values have been generated in Octave script.
    constant chrp_ampl_c : chrp_rom_t:=(
        x"F1A1",
        x"03D6",
        x"0950",
        x"EE5C",
        x"0EA5",
        x"FF51",
        x"F1A1",
        x"10FE",
        x"FD0E",
        x"F15B",
        x"0FA3",
        x"02F2",
        x"EE1B",
        x"074B",
        x"0EA5",
        x"F31F",
        x"F4D1",
        x"0EE8",
        x"0A15",
        x"F15B",
        x"F41D",
        x"0BE3",
        x"0FA3",
        x"FAD7",
        x"EDE6",
        x"F9FA",
        x"0D81",
        x"10A8",
        x"0197",
        x"F15B",
        x"EF02",
        x"FAD7",
        x"0A15",
        x"11E5",
        x"0EA5",
        x"03D6",
        x"F7E2",
        x"EFEC",
        x"EDFD",
        x"F15B",
        x"F7E2",
        x"FF51",
        x"0606",
        x"0B2F",
        x"0EA5",
        x"10A8",
        x"11A4",
        x"1203",
        x"1218",
        x"121A"
    );

    constant chrp_phs_c : chrp_rom_t:=(
        x"0B01",
        x"EE4F",
        x"0F86",
        x"FBF1",
        x"F55C",
        x"1217",
        x"F4FF",
        x"F9C3",
        x"11DC",
        x"F55C",
        x"F6E2",
        x"11DC",
        x"FD48",
        x"EF6F",
        x"0AA4",
        x"0CB8",
        x"F1C4",
        x"F5BB",
        x"0F09",
        x"0AA4",
        x"F259",
        x"F259",
        x"091E",
        x"115A",
        x"0000",
        x"EEEE",
        x"F3F2",
        x"0716",
        x"1208",
        x"0AA4",
        x"F9C3",
        x"EEA6",
        x"F0F7",
        x"FD48",
        x"0AA4",
        x"11B1",
        x"102E",
        x"0852",
        x"FE2F",
        x"F55C",
        x"EFD2",
        x"EDE9",
        x"EEEE",
        x"F1C4",
        x"F55C",
        x"F8EA",
        x"FBF1",
        x"FE2F",
        x"FF8C",
        x"0000"
    );

    -- This function rounds fix-point numbers to even(convergent rounding)
    function conv_round(a: in std_logic_vector; frc_len: in integer) return std_logic_vector;

end package radar_pkg;

package body radar_pkg is

function conv_round(a: in std_logic_vector; frc_len: in integer) return std_logic_vector is
    constant NEG_MAX    : std_logic_vector(a'length-1 downto 0) :=(a'high => '1', others => '0');
    constant POS_MAX    : std_logic_vector(a'length-1 downto 0) :=(a'high => '0', others => '1');
    variable ret_val    : std_logic_vector(a'length downto 0);
    constant frc_half   : unsigned(frc_len-1 downto 0) := shift_left(to_unsigned(1,frc_len), frc_len-1);
    variable orig_sign  : std_logic;
    variable ovf_sign   : std_logic:='0';
begin
    ovf_sign    := '0';
    orig_sign   := a(a'left);
    ret_val     := orig_sign & a;
    if frc_len > 0 then     -- fixed-point value (fractional bits present)
        if unsigned(ret_val(frc_len-1 downto 0)) = frc_half then    -- check if it is equal to 0.5!
            if ret_val(frc_len) = '1' then
                ret_val(ret_val'left downto frc_len) := std_logic_vector(unsigned(ret_val(ret_val'left downto frc_len))+1);
            end if;
        else
            ret_val := std_logic_vector(unsigned(ret_val) + frc_half);
        end if;
        if ret_val(ret_val'left-1) /= ret_val(ret_val'left) then    -- sign bit has changed
            if orig_sign = '0' then
                ret_val(a'left downto 0) := POS_MAX;
            else
                ret_val(a'left downto 0) := NEG_MAX;
            end if;
            ovf_sign    := '1';
        else
            ret_val(frc_len-1 downto 0) := (others => '0');
        end if;
    end if;
    return ovf_sign & ret_val(a'left downto 0);
end function conv_round;

end package body radar_pkg;