function [h_abs, k, w] = freq_mag_response(p, z)
    w = linspace(0,pi,320);

    for i=numel(p):-1:1
        if abs(imag(p(i))) < 1E-3
            p(i) = real(p(i));
        end
        if isreal(p(i))
            Fp_abs(i,:) = sqrt(1 - 2*p(i)*cos(w) + p(i)^2);
        elseif imag(p(i)) > 1E-3  % Skip negative part of complex pole pair
            r = abs(p(i));
            theta = angle(p(i));
            Fp_abs(i,:) = sqrt((cos(2*w) - 2*r*cos(theta)*cos(w) + r^2).^2 + (sin(2*w) - 2*r*cos(theta)*sin(w)).^2);
        end
    end
    for i=numel(z):-1:1
        if abs(imag(z(i))) < 1E-3
            z(i) = real(z(i));
        end
        if isreal(z(i))
            Fz_abs(i,:) = sqrt(1 - 2*z(i)*cos(w) + z(i)^2);
        elseif imag(z(i)) > 1E-3  % Skip negative part of complex zero pair
            r = abs(z(i));
            theta = angle(z(i));
            Fz_abs(i,:) = sqrt((cos(2*w) - 2*r*cos(theta)*cos(w) + r^2).^2 + (sin(2*w) - 2*r*cos(theta)*sin(w)).^2);
        end
    end

    h_abs = log10(prod(Fz_abs,1)./prod(Fp_abs,1))*20;
    k = max(h_abs);
    h_abs = h_abs - k;

    if nargout ~= 3
        plot(w/pi, h_abs)
        xlabel("Normalized Frequency (pi*rad/sample)")
        ylabel("Magnitude (dB)")
        title("Peak gain = " + k + "dB")
        grid on
    end
end