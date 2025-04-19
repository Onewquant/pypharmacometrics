# EBE for TDM
e = environment()

mGrad = function(func, x, nRec){
  
  # multiple gradient function in case of func returns a vector for each x (vector) point
  # Returns a nRec * length(x) matrix of gradients
  # Each row_i is the corresponding gradient of x_i
  n = length(x)
  x1 = vector(length=n)
  x2 = vector(length=n)
  ga = matrix(nrow=nRec, ncol=4)
  gr = matrix(nrow=nRec, ncol=n)
  
  for (i in 2:n) x1[i] = x2[i] = x[i]
  
  for (i in 1:n) {
    axi = abs(x[i])
    if (axi <= 1) hi = 1e-4
    else          hi = 1e-4 * axi
    
    for (k in 1:4) {
      x1[i] = x[i] - hi
      x2[i] = x[i] + hi
      ga[,k] = (func(x2) - func(x1)) / (2*hi)
      hi = hi / 2
    }
    
    ga[,1] = (ga[,2]*4 - ga[,1]) / 3
    ga[,2] = (ga[,3]*4 - ga[,2]) / 3
    ga[,3] = (ga[,4]*4 - ga[,3]) / 3
    ga[,1] = (ga[,2]*16 - ga[,1]) / 15
    ga[,2] = (ga[,3]*16 - ga[,2]) / 15
    ga[,i] = (ga[,2]*64 - ga[,1]) / 63
    x1[i] = x2[i] = x[i]
  }
  
  return(gr)
  }


objEta = function(ETAi){
  
  # External Variable : e$DATAi, e$TH, e$invOH, e$SG
  # External Function : e$PRED
  e$Fi = e$PRED(e$TH, ETAi, e$DATAi)[!is.na(e$DATAi$DV)]
  e$Ri = e$DATAi[!is.na(e$DATAi$DV), "DV"] - e$Fi
  e$Hi = cbind(e$Fi, 1)
  e$Vi = diag(e$Hi %*% e$SG %*% t(e$Hi))
  sum(log(e$Vi) + e$Ri*e$Ri/e$Vi) + t(ETAi) %*% e$invOM %*% ETAi
  }


EBE = function(PRED, DATAi, TH, OM, SG){

  # PRED = PredVanco
  # DATAi = data.frame(
  #   ID = c(1, 1, 1, 1, 1, 1, 1),
  #   TIME = c(0.00000, 11.00000, 23.00000, 35.00000, 45.83333, 46.00000, 49.33333),
  #   AMT = c(751.5, 751.5, 751.5, 751.5, 751.5, 762.0, NA),
  #   RATE = c(751.5, 751.5, 751.5, 751.5, 751.5, 762.0, NA),
  #   DV = c(NA, NA, NA, NA, 5.0, NA, 15.5),
  #   MDV = c(1, 1, 1, 1, 0, 1, 0),
  #   SEX = c(1, 1, 1, 1, 1, 1, 1),
  #   AGE = c(65, 65, 65, 65, 65, 65, 65),
  #   BWT = c(50.1, 50.1, 50.1, 50.1, 50.1, 50.8, 50.8),
  #   SCR = c(0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56),
  #   CLCR = c(93.19, 93.19, 93.19, 93.19, 93.19, 94.49, 94.49)
  # )
  # 
  # TH = c(3.8135955291021233, 39.889510090195238, 44.981835351176571, 2.0055189192561507)         # CL, V1, V2, Q?
  OM = matrix(c( 0.10855133849022583,       -1.52445093837639736E-002,  -0.19698189309256298,       0.11914555131547180,
                 -1.52445093837639736E-002, 3.01016276351715687E-003,   2.31285256338822770E-002,   -6.61597586201947800E-003,
                 -0.19698189309256298,      2.31285256338822770E-002,   1.0420667710781930,         -0.19223488058085114,
                 0.11914555131547180,       -6.61597586201947800E-003,  -0.19223488058085114,       0.416
  ), nrow = 4, byrow = TRUE)
  # 
  # SG = matrix(c(0.14019731106615912^2, 0.00000000,
  #               0.00000000,            1.8662549226759475^2), nrow=2)
  
  # e (환경)에 각 변수 저장해놓고 사용
  e$PRED = PRED
  e$DATAi = DATAi
  e$TH = TH
  e$OM = OM
  e$SG = SG
  e$nEta = nrow(OM)
  e$invOM = solve(OM)
  
  r0 = optim(rep(0, e$nEta), objEta, method="BFGS", hessian=TRUE)
# r0 = nlm(ObjEta, rep(0, e$nEta), hessian=TRUE)
  EBEi = r0$par
  COV = 2*solve(r0$hessian)
  SE = sqrt(diag(COV))
  Fi = PRED(TH, EBEi, DATAi)
  Ri = DATAi[!is.na(DATAi$DV),"DV"] - Fi[!is.na(DATAi$DV)]
  nRec = length(Fi)
  
  PREDij = function(ETA) PRED(TH, ETA, DATAi)
  gr1 = mGrad(PREDij, EBEi, nRec)   # gr is a (nRec x nEta) matrix of gradients,
  VF = diag(gr1 %*% COV %*% t(gr1)) # Transpose of gr is switched, because, gr is a matrix n~~
  SEy = sqrt(VF) # Considering only SE of EBE. SE for mean y of the same EBE individuals
  SDy = sqrt(VF + VF*SG[1,1] + Fi*Fi*SG[1,1] + SG[2,2]) # Above + Variance of epsilons. in case ~~
# Cf. If X and Y are independent, V(XY) = V(X)*V(Y) + E(X)^2*V(Y) + V(X)*E(Y)^2 (See MGB 3e.)
# for Prediction Interval considering SE of EBE and Variances of Epsilons
  
  Res = list(EBEi, SE, COV, Fi, SEy, SDy, Ri)
  names(Res) = c("EBEi", "SE", "COV", "IPRED", "SE.IPRED", "SD.IPRED", "IRES")
  return(Res)
  }


calcPI = function(PRED, DATAi, TH, SG, rEBE, npoints=500)
  {
  EBEi = rEBE$EBEi
  nEta = length(EBEi)
  COV = rEBE$COV
  DATAi2 = merge(DATAi, data.frame(TIME = seq(0, max(DATAi$TIME), length=npoints)), by="TIME")     # [유추필요] by="TIME", ~~)
  y2 = PRED(TH, EBEi, DATAi2)
  nRec2 = length(y2)
  
  # PREDij = function(ETA) PRED(TH, ETA, DATAi2)
  # gr2 = mGrad(PREDij, EBEi, nRec2)   # gr is a (nRec x nEta) matrix of gradients.
  # VF2 = diag(gr2 %*% COV %*% t(gr2)) # Transpose of gr is switched, because, gr is a matrix ~~
  # SEy2 = sqrt(VF2)
  # SDy2 = sqrt(VF2 + VF2*SG[1,1] + y2*y2*SG[1,1] + SG[2,2])
  
  
  # 예측 확인
  cat("✅ y2 range: ", range(y2, na.rm = TRUE), "\n")
  
  PREDij = function(ETA) PRED(TH, ETA, DATAi2)
  gr2 = mGrad(PREDij, EBEi, nRec2)
  
  if (any(is.na(gr2))) {
    warning("❌ Gradient (gr2) contains NA values.")
    return(NULL)
  }
  
  VF2 = diag(gr2 %*% COV %*% t(gr2))
  cat("✅ VF2 range: ", range(VF2, na.rm = TRUE), "\n")
  
  SEy2 = sqrt(VF2)
  SDy2 = sqrt(VF2 + VF2 * SG[1, 1] + y2^2 * SG[1, 1] + SG[2, 2])
  
  cat("✅ SDy2 range: ", range(SDy2, na.rm = TRUE), "\n")

  
    
  x = DATAi2$TIME
  y = DATAi2$DV
  yciLL = y2 - 1.96*SEy2
  yciUL = y2 + 1.96*SEy2
  ypiLL = y2 - 1.96*SDy2
  ypiUL = y2 + 1.96*SDy2
  
  return(data.frame(x, y, y2, yciLL, yciUL, ypiLL, ypiUL))
  }


plotPI = function(PRED, DATAi, TH, SG, rEBE, npoints=500){
  
  Res = calsPI(PRED, DATAi, TH, SG, rEBE, npoints)
  x = Res$x
  y = Res$y
  y2 = Res$y2
  yciLL = Res$yciLL
  yciUL = Res$yciUL
  ypiLL = Res$ypiLL
  ypiUL = Res$ypiUL
  dev.new()
  plot(0, 0, type="n", xlab="Time", ylab="Concentration +/- 2SD", xlim=c(min(x), max(x)), ylim=, xlim=c(min(y), max(y))) # [유추필요] xlim= ~~)
  points(x[!is.na(y)], y[!is.na(y)], pch=16)
  lines(x, y2, lty=1)
# polygon(c(x, rev(x)), c(ypiLL, rev(ypiUL)), col=#111111)
# polygon(c(x, rev(x)), c(yciLL, rev(yciUL)), col=#222222)
  lines(x, yciLL, lty=2, col="red")
  lines(x, yciLL, lty=2, col="red")
  lines(x, ypiUL, lty=2, col="blue")
  lines(x, ypiUL, lty=2, col="blue")
  abline(h=c(5,15,25,35, lty=2))
  return(Res)
  }



calcTDM = function(PRED, DATAi, TH, SG, rEBE, TIME, AMT, RATE, II, ADDL, npoints=500)
{
  DATAi = addDATAi(DATAi, TIME, AMT, RATE, II, ADDL)
  rTab = calcPI(PRED, DATAi, TH, SG, rEBE, npoints)
  return(rTab)
}


plotTDM = function(PRED, DATAi, TH, SG, rEBE, TIME, AMT, RATE, II, ADDL, npoints=500)
{
  DATAi = addDATAi(DATAi, TIME, AMT, RATE, II, ADDL)
  rTab = plotPI(PRED, DATAi, TH, SG, rEBE, npoints)
  return(rTab)
}


convDT = function(DATAi)
{
  DATAi[,"SDT"] = strftime(paste(DATAi[,"DATE"], DATAi[,"TIME"]))
  FDT = DATAi[1, "SDT"]
  DATAi[,"TIME"] = as.numeric(difftime(DATAi[,"SDT"], FDT, unit="hours"))
  return(DATAi[,setdiff(colnames(DATAi),"DATE")])
}


addDATAi = function(DATAi, TIME, AMT, RATE, II, ADDL)
{
  lRow = DATAi[nrow(DATAi),]

  nADD = ADDL + 1
  for (i in 1:nADD) {
    aRow = lRow
    aRow[,"TIME"] = TIME + (i - 1)*II
    aRow[,"AMT"]  = AMT
    aRow[,"RATE"] = RATE
    aRow[,"DV"]   = NA
    aRow[,"MDV"]  = 1
    DATAi = rbind(DATAi, aRow)
  }
  return(DATAi)
}

# addDATAi = function(DATAi, TIME, AMT, RATE, II, ADDL) {
#   lRow = DATAi[nrow(DATAi), ]
#   
#   if (is.null(ADDL) || is.na(ADDL || ADDL==0)) {
#     nADD = 0
#   } else {
#     nADD = ADDL
#   }
#   
#   for (i in 0:nADD) {
#     aRow = lRow
#     aRow["TIME"] = TIME + i * II
#     aRow["AMT"]  = AMT
#     aRow["RATE"] = RATE
#     aRow["DV"]   = NA
#     aRow["MDV"]  = 1
#     DATAi = rbind(DATAi, aRow)
#   }
#   
#   return(DATAi)
# }


expandDATA = function(DATAo){
  
  eDATAi =data.frame()
  Added = vector()
  for (i in 1:nrow(DATAo)) {
    if (!is.na(DATAo[i,"ADDL"])) {
      eDATAi = rbind(eDATAi, DATAo[i,])
      cTIME  = DATAo[i,"TIME"]
      cII    = DATAo[i,"II"]
      Added = c(Added, FALSE)
      nADD = DATAo[i, "ADDL"]
      cRow = DATAo[i,]
      for (j in 1:nADD) {
        cRow[,"TIME"] = cTIME + j*cII
        eDATAi = rbind(eDATAi, cRow)
        Added = c(Added, TRUE)
      }
    } else{
      eDATAi = rbind(eDATAi, DATAo[i,])
      Added = c(Added, FALSE)
    }
  }
  eDATAi[,"II"] = NA
  eDATAi[,"ADDL"] = NA
  
  AddOrder = Added[order(eDATAi$TIME)]
  Added <<- AddOrder
  eDATAi = eDATAi[order(eDATAi$TIME),]
  
  CovCol = setdiff(colnames(eDATAi), c("ID","TIME","AMT","RATE","II","ADDL","DV","MDV")) # [유추필요] ,"DV","MDV",~~))
  for (i in 2:nrow(eDATAi)) {
    if (AddOrder[i] == TRUE) eDATAi[i, CovCol] = eDATAi[i -1, CovCol]
  }
  
  return(eDATAi)
}


PredVanco = function(TH, ETA, DATAi){

  V1 = TH[2]*exp(ETA[2])
  V2 = TH[3]*exp(ETA[3])
  Q = TH[4]*exp(ETA[4])
  K12 = Q/V1
  K21 = Q/V2

  IPRE = 0
  inf = FALSE
  Conc = data.frame(C1=0, C2=0)
  pCLCR = min(DATAi[1,"CLCR"], 150)
  for (i in 2:nrow(DATAi)) {
    CLCR = min(DATAi[i,"CLCR"], 150)
    if (is.na(CLCR)) {
      CLCR = pCLCR
    } else {
      pCLCR = CLCR
    }
    TVCL = TH[1]*CLCR/100
    CL = TVCL*exp(ETA[1])
    K10 = CL/V1
    AlpBe = K10 + K12 + K21
    AlmBe = K10*K21
    Det4 = sqrt(AlpBe*AlpBe - 4*AlmBe)
    Alpha = (AlpBe + Det4)/2
    Beta = (AlpBe - Det4)/2
    Divisor = V1*(Alpha - Beta)
    cTime = DATAi[i, "TIME"]
    dTime = cTime - DATAi[i - 1, "TIME"]

    if (!is.na(DATAi[i-1, "AMT"]) & DATAi[i-1,"AMT"] > 0) { # inf should not be TRUE
      pTime = DATAi[i-1,"TIME"]
      pAMT = DATAi[i-1,"AMT"]
      pRate = DATAi[i-1,"RATE"]
      pDur = pAMT/pRate
      if (dTime <= pDur) { # Infusion times should not overlab
        infC1 = pRate*(Alpha - K21)/Divisor*(1 - exp(-Alpha*dTime))/Alpha
        infC2 = pRate*(K21 - Beta)/Divisor*(1 - exp(-Beta*dTime))/Beta
        C1 = Conc[i - 1, "C1"]*exp(-Alpha*dTime) + infC1
        C2 = Conc[i - 1, "C2"]*exp(-Beta*dTime) + infC2
        Conc = rbind(Conc, c(C1, C2))
        inf = TRUE
      } else {
        eC1 = pRate*(Alpha - K21)/Divisor*(1 - exp(-Alpha*pDur))/Alpha
        eC2 = pRate*(K21 - Beta)/Divisor*(1 - exp(-Beta*pDur))/Beta
        pC1 = Conc[i-1,"C1"]*exp(-Alpha*pDur)
        pC2 = Conc[i-1,"C2"]*exp(-Beta*pDur)
        C1 = (eC1 + pC1)*exp(-Alpha*(dTime - pDur))
        C2 = (eC2 + pC2)*exp(-Beta*(dTime - pDur))
        Conc = rbind(Conc, c(C1, C2))
        inf = FALSE
        }
      } else if (inf == TRUE) { # Infusion has not finished. Use last infusion information
      if (cTime <= pTime + pDur) {
        infC1 = pRate*(Alpha - K21)/Divisor*(1 - exp(Alpha*dTime))/Alpha
        infC2 = pRate*(K21 - Beta)/Divisor*(1 - exp(Beta*dTime))/Beta
        C1 = Conc[i-1, "C1"]*exp(-Alpha*dTime) + infC1
        C2 = Conc[i-1, "C2"]*exp(-Beta*dTime) + infC2
        Conc = rbind(Conc, c(C1, C2))
        inf = TRUE # no change of pTime, pAMT, pRate, pDur
      } else {
        rDur = pTime + pDur - DATAi[i - 1, "TIME"]
        eC1 = pRate*(Alpha - K21)/Divisor*(1-exp(-Alpha*rDur))/Alpha
        eC2 = pRate*(K21 - Beta)/Divisor*(1-exp(-Beta*rDur))/Beta
        pC1 = Conc[i - 1, "C1"]*exp(-Alpha*rDur)
        pC2 = Conc[i - 1, "C2"]*exp(-Beta*rDur)
        C1 = (eC1 + pC1) * exp(-Alpha*(cTime - pTime - spDur))
        C2 = (eC2 + pC2) * exp(-Beta*(cTime - pTime - pDur))
        Conc = rbind(Conc, c(C1, C2))
        inf = FALSE
      }

      } else {
        C1 = Conc[i - 1, "C1"]*exp(-Alpha*dTime)
        C2 = Conc[i - 1, "C2"]*exp(-Beta*dTime)
        Conc = rbind(Conc, c(C1, C2))
        inf = FALSE
      }
      IPRE[i] = C1 + C2
    }
  return(IPRE)
  }

# PredVanco = function(TH, ETA, DATAi) {
#   V1 = TH[2] * exp(ETA[2])
#   V2 = TH[3] * exp(ETA[3])
#   Q  = TH[4] * exp(ETA[4])
#   K12 = Q / V1
#   K21 = Q / V2
#   
#   IPRE = rep(NA, nrow(DATAi))  # 미리 NA로 초기화
#   inf = FALSE
#   Conc = data.frame(C1 = 0, C2 = 0)
#   
#   pCLCR = min(DATAi[1, "CLCR"], 150)
#   
#   for (i in 2:nrow(DATAi)) {
#     CLCR = min(DATAi[i, "CLCR"], 150)
#     if (is.na(CLCR)) {
#       CLCR = pCLCR
#     } else {
#       pCLCR = CLCR
#     }
#     
#     TVCL = TH[1] * CLCR / 100
#     CL = TVCL * exp(ETA[1])
#     K10 = CL / V1
#     
#     AlpBe = K10 + K12 + K21
#     AlmBe = K10 * K21
#     Det4_expr = AlpBe^2 - 4 * AlmBe
#     
#     # 안정성 체크
#     if (!is.finite(Det4_expr) || Det4_expr < 0 || !is.finite(CL) || !is.finite(K10)) {
#       warning(paste("❌ Invalid root or clearance at row", i, "| CLCR:", CLCR, "| ETA:", paste(ETA, collapse = ",")))
#       return(rep(NA, nrow(DATAi)))
#     }
#     
#     Det4 = sqrt(Det4_expr)
#     Alpha = (AlpBe + Det4) / 2
#     Beta  = (AlpBe - Det4) / 2
#     Divisor = V1 * (Alpha - Beta)
#     
#     if (!is.finite(Divisor) || Divisor == 0) {
#       warning(paste("❌ Invalid divisor at row", i, "| Divisor:", Divisor))
#       return(rep(NA, nrow(DATAi)))
#     }
#     
#     cTime = DATAi[i, "TIME"]
#     dTime = cTime - DATAi[i - 1, "TIME"]
#     
#     if (!is.na(DATAi[i - 1, "AMT"]) && DATAi[i - 1, "AMT"] > 0) {
#       pTime = DATAi[i - 1, "TIME"]
#       pAMT  = DATAi[i - 1, "AMT"]
#       pRate = DATAi[i - 1, "RATE"]
#       pDur  = pAMT / pRate
#       
#       if (dTime <= pDur) {
#         infC1 = pRate * (Alpha - K21) / Divisor * (1 - exp(-Alpha * dTime)) / Alpha
#         infC2 = pRate * (K21 - Beta) / Divisor * (1 - exp(-Beta * dTime)) / Beta
#         C1 = Conc[i - 1, "C1"] * exp(-Alpha * dTime) + infC1
#         C2 = Conc[i - 1, "C2"] * exp(-Beta * dTime) + infC2
#         Conc = rbind(Conc, c(C1, C2))
#         inf = TRUE
#       } else {
#         eC1 = pRate * (Alpha - K21) / Divisor * (1 - exp(-Alpha * pDur)) / Alpha
#         eC2 = pRate * (K21 - Beta) / Divisor * (1 - exp(-Beta * pDur)) / Beta
#         pC1 = Conc[i - 1, "C1"] * exp(-Alpha * pDur)
#         pC2 = Conc[i - 1, "C2"] * exp(-Beta * pDur)
#         C1 = (eC1 + pC1) * exp(-Alpha * (dTime - pDur))
#         C2 = (eC2 + pC2) * exp(-Beta * (dTime - pDur))
#         Conc = rbind(Conc, c(C1, C2))
#         inf = FALSE
#       }
#       
#     } else if (inf == TRUE) {
#       if (cTime <= pTime + pDur) {
#         infC1 = pRate * (Alpha - K21) / Divisor * (1 - exp(-Alpha * dTime)) / Alpha
#         infC2 = pRate * (K21 - Beta) / Divisor * (1 - exp(-Beta * dTime)) / Beta
#         C1 = Conc[i - 1, "C1"] * exp(-Alpha * dTime) + infC1
#         C2 = Conc[i - 1, "C2"] * exp(-Beta * dTime) + infC2
#         Conc = rbind(Conc, c(C1, C2))
#         inf = TRUE
#       } else {
#         rDur = pTime + pDur - DATAi[i - 1, "TIME"]
#         eC1 = pRate * (Alpha - K21) / Divisor * (1 - exp(-Alpha * rDur)) / Alpha
#         eC2 = pRate * (K21 - Beta) / Divisor * (1 - exp(-Beta * rDur)) / Beta
#         pC1 = Conc[i - 1, "C1"] * exp(-Alpha * rDur)
#         pC2 = Conc[i - 1, "C2"] * exp(-Beta * rDur)
#         C1 = (eC1 + pC1) * exp(-Alpha * (cTime - pTime - pDur))
#         C2 = (eC2 + pC2) * exp(-Beta * (cTime - pTime - pDur))
#         Conc = rbind(Conc, c(C1, C2))
#         inf = FALSE
#       }
#       
#     } else {
#       C1 = Conc[i - 1, "C1"] * exp(-Alpha * dTime)
#       C2 = Conc[i - 1, "C2"] * exp(-Beta * dTime)
#       Conc = rbind(Conc, c(C1, C2))
#       inf = FALSE
#     }
#     
#     if (!is.finite(C1 + C2)) {
#       warning(paste("❌ Non-finite prediction at row", i, "→ C1+C2 not finite"))
#       return(rep(NA, nrow(DATAi)))
#     }
#     
#     IPRE[i] = C1 + C2
#   }
#   
#   return(IPRE)
# }