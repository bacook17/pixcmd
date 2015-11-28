FUNCTION ADD_PHOT_ERR, str, i1,i2, p1,p2

  err = 10^(p1[0]+p1[1]*str.(i1))
  str.(i1) = str.(i1) + randomn(seed,n_elements(str))*err

  err = 10^(p2[0]+p2[1]*str.(i2))
  str.(i2) = str.(i2) + randomn(seed,n_elements(str))*err

  RETURN, str

END
