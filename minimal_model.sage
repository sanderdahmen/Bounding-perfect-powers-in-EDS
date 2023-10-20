def minimal_maker(E, u, a1, a2, a3):
    s = (u*a1 - E.a1())/2
    r = (u^2*a2 - E.a2() + s*E.a1() + s^2)/3
    t = (u^3*a3 - E.a3() - r*E.a1())/2
    a4 = (E.a4() - s*E.a3() + 2*r*E.a2() - (t + r*s)*E.a1() + 3*r^2 - 2*s*t)/(u^4)
    a6 = (E.a6() + r*E.a4() + r^2*E.a2() + r^3 - t*E.a3() - t^2 - r*t*E.a1())/(u^6)
    return EllipticCurve([a1, a2, a3, a4, a6])

