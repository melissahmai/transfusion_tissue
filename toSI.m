function converted = toSI(toconvert,h,eta,D,cm)
    converted = toconvert;
    converted.C = toconvert.C*cm;
    converted.Cm = toconvert.Cm*cm;
    converted.Cp = toconvert.Cp*cm;
    converted.P = toconvert.P*eta*D/(h^2);
    converted.Pm = toconvert.Pm*eta*D/(h^2);
    converted.Psi = toconvert.Psi*eta*D/(h^2);
    converted.Uij = toconvert.Uij*D/h;
    converted.Jij = toconvert.Jij*D*cm/h;
    converted.Uc = toconvert.Uc*D*cm/h;
    converted.Dc = toconvert.Dc*D*cm/h;
    converted.Uib = toconvert.Uib*D/h;
    converted.Jib = toconvert.Jib*D*cm/h;
    converted.RT = toconvert.RT*eta*D/(cm*h^2);
end