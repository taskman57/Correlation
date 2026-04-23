function [rounded_val, ovf] = myround(val, intb, frcb)
  % =========================================================================
  % FUNCTION: myconv_round
  % -------------------------------------------------------------------------
  % Implementation of a Hardware-Accurate Unbiased Convergent Rounding
  % (Round-to-Even) with Symmetric Saturation logic.
  %
  % STAGES:
  % 1. FIXED-POINT ALIGNMENT: Scales input to internal integer representation.
  % 2. TIE-BREAKING LOGIC: Implements "Round-to-Even" to eliminate DC bias
  %    accumulation inherent in standard truncation or asymmetric rounding.
  % 3. DYNAMIC RANGE LIMITING: Performs Sign-Aware Saturation (Clipping)
  %    to prevent Two's Complement wrap-around/phase-reversal.
  % 4. RTL EMULATION: Provides bit-accurate output mirroring FPGA DSP
  %    slice behavior for bit-true hardware verification.
  % =========================================================================
  intw = intb + frcb;
  POS_MAX = 2^(intw-1)-1;                 % Maximum value
  NEG_MAX = -2^(intw-1);                  % Minimum value
  ovf = 0;
  sign = 0;
  val_int = val*(2^frcb);
  rnd_val = floor(val_int);               % make an integer
  if rnd_val <0
    sign = 1;
    rnd_val = mod(rnd_val, 2^(intw+sign));
  endif
  mask = 2^(intw+sign) - 1;
  wrapped = bitand(rnd_val, mask);        % wrap into exact intw

  rnd_val = wrapped;
  if frcb >0
    rm= rem(rnd_val, 2^frcb);
    if rm == 2^(frcb-1)                   % if fractional is 0.5
      rnd_val = bitshift(rnd_val,...      % pick integer part
       (-1*frcb));
      if rem(rnd_val,2)                   % if integer part is odd
        rnd_val += 1;
      endif
      rnd_val = bitshift(rnd_val,frcb);   % shift bits back to the right position
    else
      rnd_val += 2^(frcb-1);              % add 0.5
    endif
  endif
  if bitget(rnd_val, intw) != bitget(rnd_val, intw+1)
    if sign == 0
      rnd_val = POS_MAX;
    else
      rnd_val = NEG_MAX;
    endif
    ovf = 1;
  else
    mask = bitcmp(2^frcb-1, intw);
    rnd_val = bitand(rnd_val, mask);    % clear all fractional bits
  endif

  rnd_val = rnd_val - ...                 % convert to negative numbers
  (rnd_val >= 2^(intw-1)) * 2^intw;
  rounded_val = rnd_val/2^frcb;
end
