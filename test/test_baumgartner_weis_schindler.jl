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

test("BWS statistic") do

  test("a small hand-calculated BWS case with 2 elements in each sample") do

    x = [1.0, 2.0]
    n = 2
    y = [3.0, 4.0]
    m = 2

    r = [1, 2] # ranks for values in x
    s = [3, 4] # ranks for values in y

    mn = m+n
    mp1 = m+1
    np1 = n+1

    b_x_part_for_i_1 = (r[1] - 1*mn/n)^2 / ( (1/np1) * (1 - 1/np1) * (m * mn / n) )
    b_x_part_for_i_2 = (r[2] - 2*mn/n)^2 / ( (2/np1) * (1 - 2/np1) * (m * mn / n) )
    b_x = (1 / n) * ( b_x_part_for_i_1 + b_x_part_for_i_2 )

    b_y_part_for_i_1 = (s[1] - 1*mn/m)^2 / ( (1/mp1) * (1 - 1/mp1) * (n * mn / m) )
    b_y_part_for_i_2 = (s[2] - 2*mn/m)^2 / ( (2/mp1) * (1 - 2/mp1) * (n * mn / m) )
    b_y = (1 / m) * ( b_y_part_for_i_1 + b_y_part_for_i_2 )

    b = (b_x + b_y) / 2

    @t FeldtLib.baumgartner_weis_schindler_statistic(x, y) == b

  end

  test("a small hand-calculated BWS case with 3+2 elemnts") do

    x = [1.0, 2.0, 5.0]
    n = 3
    y = [3.0, 4.0]
    m = 2

    r = [1, 2, 5] # ranks for values in x
    s = [3, 4] # ranks for values in y

    mn = m+n
    mp1 = m+1
    np1 = n+1

    bx1 = (r[1] - 1*mn/n)^2 / ( (1/np1) * (1 - 1/np1) * (m * mn / n) )
    bx2 = (r[2] - 2*mn/n)^2 / ( (2/np1) * (1 - 2/np1) * (m * mn / n) )
    bx3 = (r[3] - 3*mn/n)^2 / ( (3/np1) * (1 - 3/np1) * (m * mn / n) )
    bx = (1 / n) * ( bx1 + bx2 + bx3)

    by1 = (s[1] - 1*mn/m)^2 / ( (1/mp1) * (1 - 1/mp1) * (n * mn / m) )
    by2 = (s[2] - 2*mn/m)^2 / ( (2/mp1) * (1 - 2/mp1) * (n * mn / m) )
    by = (1 / m) * ( by1 + by2 )

    b = (bx + by) / 2

    @t FeldtLib.baumgartner_weis_schindler_statistic(x, y) == b

  end

end