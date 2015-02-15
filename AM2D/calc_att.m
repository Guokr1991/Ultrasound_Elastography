function xx = calc_att (a_t, f_0, f_s)

% a_t : frequency dependent attenuation coefficient, in dB/cm/MHz
% f_s : ultrasound center freq. in MHz
% f_s : sampling freq. in MHz
% xx : x in the paper Equation 18
xx = exp(1540e2*a_t*f_0*log(10)/20e6/f_s); % to compensate for attenuation