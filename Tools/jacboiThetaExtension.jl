using EllipticFunctions

function th1(z, tau)
    return jtheta1(z, exp(im * pi * tau))
end
function th2(z, tau)
    return jtheta2(z, exp(im * pi * tau))
end

function th3(z, tau)
    return jtheta3(z, exp(im * pi * tau))
end

function th4(z, tau)
    return jtheta4(z, exp(im * pi * tau))
end
function th1prime(z, tau)
    return jtheta1dash(z, exp(im * pi * tau))
end
function th2prime(z, tau)
    return jtheta1dash(z + 0.5 * pi, exp(im * pi * tau))
end
function th4prime(z, tau)
    q = exp(im * pi * tau)
    coeff = exp(-im) * q^(0.25)
    
    return coeff * (jtheta1(z + 0.5 * im * log(q), q) + im * jtheta1dash(z + 0.5 * im * log(q), q))
end