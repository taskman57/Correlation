function [rounded_val, ovf] = my_comp_round(val, intb, frcb)
  if ~isreal(val)
      % Process Real and Imaginary parts separately
      [r_val, r_ovf] = myround(real(val), intb, frcb);
      [i_val, i_ovf] = myround(imag(val), intb, frcb);

      rounded_val = r_val + 1i * i_val;
      ovf = r_ovf | i_ovf;
  else
      % Fallback for real numbers
      [rounded_val, ovf] = myround(val, intb, frcb);
  end
end
