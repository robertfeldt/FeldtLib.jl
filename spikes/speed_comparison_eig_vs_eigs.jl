function gen_symmetric_matrix(n)
  a = randn(n,n)
  (a + a')/2
end

for n in [10, 100, 500, 1000, 2000]
  a = gen_symmetric_matrix(n)

  tic()
    de, ve = eig(a)
  teig = toq()

  k = 5
  while k <= n/2
    tic()
      des, ves = eigs(a; nev = k)
    teigs = toq()

    #dnorm = norm( de-des[1:k] )
    #vnorm = norm( ve-ves )

    println("$(n): k = $(k), eigd norm = , eigv norm = , eigs time factor = $(teig/teigs), eig time = $(teig), eigs time = $(teigs)")

    k *= 10
  end
end