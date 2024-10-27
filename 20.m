clear

p = [0.9738*exp(pi/10*1j)];
z = [exp(pi/10*1j)];

omega = linspace(0,pi,320);

for i=numel(p):-1:1
    if isreal(p(i))
        Fp_abs(i,:) = sqrt(1 - 2*p(i)*cos(omega) + p(i)^2);
    else
        r = abs(p(i));
        theta = angle(p(i));
        Fp_abs(i,:) = sqrt((cos(2*omega) - 2*r*cos(theta)*cos(omega) + r^2).^2 + (sin(2*omega) - 2*r*cos(theta)*sin(omega)).^2);
    end
end
for i=numel(z):-1:1
    if isreal(z(i))
        Fz_abs(i,:) = sqrt(1 - 2*z(i)*cos(omega) + z(i)^2);
    else
        r = abs(z(i));
        theta = angle(z(i));
        Fz_abs(i,:) = sqrt((cos(2*omega) - 2*r*cos(theta)*cos(omega) + r^2).^2 + (sin(2*omega) - 2*r*cos(theta)*sin(omega)).^2);
    end
end

H_abs = log10(prod(Fz_abs,1)./prod(Fp_abs,1))*20;
peak_gain = max(H_abs);
H_abs = H_abs - peak_gain;

plot(omega/pi, H_abs)
xlabel("Normalized Frequency (pi*rad/sample)")
ylabel("Magnitude (dB)")
title("Peak gain = " + peak_gain + "dB")
grid on