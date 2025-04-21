# R to python project
e = environment()

raw_data = read.csv("C:/Users/ilma0/PycharmProjects/pypharmacometrics/Projects/TDM_Engine/sample_vanco_loading_maintenance.csv", na.strings = c("", ".", "NA"), as.is = TRUE)


TH = c(3.8135955291021233, 39.889510090195238, 44.981835351176571, 2.0055189192561507)         # CL, V1, V2, Q?
OM = matrix(c( 0.10855133849022583,       -1.52445093837639736E-002,  -0.19698189309256298,       0.11914555131547180,
               -1.52445093837639736E-002, 3.01016276351715687E-003,   2.31285256338822770E-002,   -6.61597586201947800E-003,
               -0.19698189309256298,      2.31285256338822770E-002,   1.0420667710781930,         -0.19223488058085114,
               0.11914555131547180,       -6.61597586201947800E-003,  -0.19223488058085114,       0.416
), nrow = 4, byrow = TRUE)

SG = matrix(c(0.14019731106615912^2, 0.00000000,
              0.00000000,            1.8662549226759475^2), nrow=2)



convDT = function(DATAi)
{
  DATAi[,"SDT"] = strftime(paste(DATAi[,"DATE"], DATAi[,"TIME"]))
  FDT = DATAi[1, "SDT"]
  DATAi[,"TIME"] = as.numeric(difftime(DATAi[,"SDT"], FDT, unit="hours"))
  return(DATAi[,setdiff(colnames(DATAi),"DATE")])
}


data_prepped <- convDT(raw_data)           # TIME 컬럼 붙여줌 (numerical time으로)
# 데이터에 ADDL 컬럼이 있고, 값이 실제로 > 0인 경우만 확장 적용
if ("ADDL" %in% names(data_prepped) && any(!is.na(data_prepped$ADDL) & data_prepped$ADDL > 0)) {
  data_prepped <- expandDATA(data_prepped)
}

DATAi = (data_prepped[, c("ID","TIME","AMT","RATE","DV","MDV","SEX","AGE","BWT","SCR","CLCR")])

# rEBE = EBE(PredVanco, DATAi, TH, OM, SG) #######################


mGrad = function(func, x, nRec){
  
  # multiple gradient function in case of func returns a vector for each x (vector) point
  # Returns a nRec * length(x) matrix of gradients
  # Each row_i is the corresponding gradient of x_i
  
  
  # func = PREDij
  # x = EBEi
  # nRec = length(Fi)
  
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
    gr[,i] = (ga[,2]*64 - ga[,1]) / 63
    x1[i] = x2[i] = x[i]
  }
  
  return(gr)
}


EBE = function(PRED, DATAi, TH, OM, SG){
  
  # PRED = PredVanco

  
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
  Fi = PRED(TH, EBEi, DATAi)      # IPRED
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


objEta = function(ETAi){
  
  # External Variable : e$DATAi, e$TH, e$invOH, e$SG
  # External Function : e$PRED
  e$Fi = e$PRED(e$TH, ETAi, e$DATAi)[!is.na(e$DATAi$DV)]
  e$Ri = e$DATAi[!is.na(e$DATAi$DV), "DV"] - e$Fi
  e$Hi = cbind(e$Fi, 1)
  e$Vi = diag(e$Hi %*% e$SG %*% t(e$Hi))
  sum(log(e$Vi) + e$Ri*e$Ri/e$Vi) + t(ETAi) %*% e$invOM %*% ETAi
}

PredVanco = function(TH, ETA, DATAi){
  
  # TH
  # ETA = EBEi
  # DATAi = DATAi2
  
  
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
        infC1 = pRate*(Alpha - K21)/Divisor*(1 - exp(-Alpha*dTime))/Alpha
        infC2 = pRate*(K21 - Beta)/Divisor*(1 - exp(-Beta*dTime))/Beta
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
        C1 = (eC1 + pC1) * exp(-Alpha*(cTime - pTime - pDur))
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




