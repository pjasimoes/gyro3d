pro freefree, temperature, np, freq, jff, kff

    ;; free-free coefficient by Dulk (1985).
    ;; assuming fully ionized hidrogen isothermal plasma.
    kb = 1.38000d-016           ;; Boltzmann constant
    c=2.9979246d+10             ;; speed of light

    t1 = (TEMPERATURE LT 2e5) ? 18.2 + alog(TEMPERATURE^(1.5)) - alog(freq) : $
         24.5 + alog(TEMPERATURE) - alog(freq)

    kff = 9.78e-3 * np^2 / freq^2 / TEMPERATURE^(1.5) * t1

    jff = kff * kb * TEMPERATURE * FREQ^2 / C^2

END 
