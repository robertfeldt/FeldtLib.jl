test("Ranks and permutations") do

  test("ranks") do

    @t FeldtLib.ranks([1]) == [1]
  
    @t FeldtLib.ranks([1, 2]) == [1, 2]
    @t FeldtLib.ranks([2, 1]) == [2, 1]
  
    @t FeldtLib.ranks([3, 2, 1]) == [3, 2, 1]
    @t FeldtLib.ranks([3, 1, 2]) == [3, 1, 2]
    @t FeldtLib.ranks([2, 1, 3]) == [2, 1, 3]
    @t FeldtLib.ranks([1, 2, 3]) == [1, 2, 3]
    @t FeldtLib.ranks([1, 3, 2]) == [1, 3, 2]
  
    @t FeldtLib.ranks([4.0, 2.5, 1.2, 5.6]) == [3, 2, 1, 4]
  
    @t FeldtLib.ranks([7.0, 3.5, 2.2, 10.6, 456.0]) == [3, 2, 1, 4, 5]

  end

  test("factorial") do
    @t factorial(1) == 1
    @t factorial(2) == 2
    @t factorial(3) == 6
    @t factorial(4) == 24
    @t factorial(5) == 120
  end

  test("combinations") do
    @t combinations(7, 14) == 3432
  end
end