FUNCTION ADD_PHOT_ERR, str,i1,i2,dm,zpt,exptime

  ;compute total counts
  ct1 = 10^(-2./5*(str.(i1)+dm-zpt[0]))*exptime[0]
  ct2 = 10^(-2./5*(str.(i2)+dm-zpt[1]))*exptime[1]

  cti1 = ct1*0.0
  cti2 = ct2*0.0

  ;draw from a Poisson distribution
  FOR i=0,n_elements(ct1)-1 DO BEGIN
     cti1[i] = poidev(ct1[i])
     cti2[i] = poidev(ct2[i])
  ENDFOR

  ;convert back to mags
  str.(i1) = -2.5*alog10(cti1/exptime[0])+zpt[0]-dm
  str.(i2) = -2.5*alog10(cti2/exptime[1])+zpt[1]-dm

  RETURN, str

END
