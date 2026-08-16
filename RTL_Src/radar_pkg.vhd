library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package radar_pkg is

    type chrp_rom_t is array(0 to 49) of signed(15 downto 0);

    -- These Inphase & quadrature values have been generated in Octave script.
    constant chrp_ampl_c : chrp_rom_t:=(
        x"F14A",
        x"049D",
        x"11E1",
        x"0A8C",
        x"F865",
        x"EE0E",
        x"F3B1",
        x"02D5",
        x"0F58",
        x"1190",
        x"09CC",
        x"FD9E",
        x"F333",
        x"EE33",
        x"EF21",
        x"F460",
        x"FB9C",
        x"02D5",
        x"08D2",
        x"0D1E",
        x"0FCF",
        x"1140",
        x"11E1",
        x"1213",
        x"121A"
    );
    constant chrp_phs_c : chrp_rom_t:=(
        x"F574",
        x"EE7F",
        x"FD2B",
        x"0EB6",
        x"106D",
        x"0262",
        x"F2BA",
        x"EE1F",
        x"F666",
        x"0464",
        x"0F39",
        x"11F2",
        x"0CCD",
        x"0348",
        x"F971",
        x"F220",
        x"EE70",
        x"EE1F",
        x"F031",
        x"F387",
        x"F72E",
        x"FA84",
        x"FD2B",
        x"FEFA",
        x"FFE3"
    );

    -- Gain of TX/RX antenna and Calculated System Parameters (R=7000m):
    -- Antenna gain for both TX & RX: 20 dBi
    -- - Pt: 1000.00 W, MF Gain: 16.99 dB
    -- - Single-Pulse SNR (Theoretical): -20.25 dB
    -- - S_watts: 3.78e-16 W, N_thermal: 4.00e-14 W
    -- Pulse 1: I=-1711, Q=-5994 | Octave Phase: -1.8489 rad | CORDIC Expected (Fixed): -4821
    -- Pulse 2: I=-2680, Q=-4921 | Octave Phase: -2.0695 rad | CORDIC Expected (Fixed): -5396
    -- Pulse 3: I=-829, Q=-7709 | Octave Phase: -1.6779 rad | CORDIC Expected (Fixed): -4375
    -- Pulse 4: I=-2485, Q=-3332 | Octave Phase: -2.2116 rad | CORDIC Expected (Fixed): -5767
    -- Pulse 5: I=-1775, Q=2701 | Octave Phase: 2.1522 rad | CORDIC Expected (Fixed): 5612
    -- Pulse 6: I=-4592, Q=1950 | Octave Phase: 2.7400 rad | CORDIC Expected (Fixed): 7145
    -- Pulse 7: I=-3360, Q=-7064 | Octave Phase: -2.0148 rad | CORDIC Expected (Fixed): -5254
    -- Pulse 8: I=-6348, Q=-8530 | Octave Phase: -2.2106 rad | CORDIC Expected (Fixed): -5764
    -- Pulse 9: I=-2013, Q=-8433 | Octave Phase: -1.8051 rad | CORDIC Expected (Fixed): -4707
    -- Pulse 10: I=-3573, Q=-4566 | Octave Phase: -2.2348 rad | CORDIC Expected (Fixed): -5827
    -- Pulse 11: I=-275, Q=-4227 | Octave Phase: -1.6358 rad | CORDIC Expected (Fixed): -4265
    type pulse_IQ_t is array(0 to 11) of signed(15 downto 0);
    signal inph_c : pulse_IQ_t := (
        to_signed(-1711,16),
        to_signed(-2680,16),
        to_signed(-829,16),
        to_signed(-2485,16),
        to_signed(-1775,16),
        to_signed(-4592,16),
        to_signed(-3360,16),
        to_signed(-6348,16),
        to_signed(-2013,16),
        to_signed(-3573,16),
        to_signed(-275,16),
        to_signed(-11357,16)
    );

    signal quadr_c : pulse_IQ_t := (
        to_signed(-5994,16),
        to_signed(-4921,16),
        to_signed(-7709,16),
        to_signed(-3332,16),
        to_signed(2701,16),
        to_signed(1950,16),
        to_signed(-7064,16),
        to_signed(-8530,16),
        to_signed(-8433,16),
        to_signed(-4566,16),
        to_signed(-4227,16),
        to_signed(-10190,16)
    );

    -- This function rounds fix-point numbers to even(convergent rounding)
    function conv_round(a: in std_logic_vector; frc_len: in integer) return std_logic_vector;

end package radar_pkg;

package body radar_pkg is

function conv_round(a: in std_logic_vector; frc_len: in integer) return std_logic_vector is
    constant NEG_MAX    : std_logic_vector(a'length-1 downto 0) :=(a'high => '1', others => '0');
    constant POS_MAX    : std_logic_vector(a'length-1 downto 0) :=(a'high => '0', others => '1');
    variable ret_val    : std_logic_vector(a'length downto 0);
    variable orig_sign  : std_logic;
    variable ovf_sign   : std_logic:='0';
begin
    ovf_sign    := '0';
    orig_sign   := a(a'left);
    ret_val     := orig_sign & a;
    if frc_len > 0 then     -- fixed-point value (fractional bits present)
        if unsigned(ret_val(frc_len-1 downto 0)) = shift_left(to_unsigned(1,frc_len), frc_len-1) then    -- check if it is equal to 0.5!
            if ret_val(frc_len) = '1' then
                ret_val(ret_val'left downto frc_len) := std_logic_vector(unsigned(ret_val(ret_val'left downto frc_len))+1);
            end if;
        else
            ret_val := std_logic_vector(unsigned(ret_val) + shift_left(to_unsigned(1,frc_len), frc_len-1));
        end if;
    end if;
    if ret_val(ret_val'left-1) /= ret_val(ret_val'left) then    -- sign bit has changed
        if orig_sign = '0' then
            ret_val(a'left downto 0) := POS_MAX;
        else
            ret_val(a'left downto 0) := NEG_MAX;
        end if;
        ovf_sign    := '1';
    else
        if frc_len > 0 then
            ret_val(frc_len-1 downto 0) := (others => '0');
        end if;
    end if;
    return ovf_sign & ret_val(a'left downto 0);
end function conv_round;

end package body radar_pkg;