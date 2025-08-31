import numpy as np
import pandas as pd

class vanco_adult():
    def __init__(self):
        # Typical parameters (THETA): CL, V1, V2, Q
        self.TH = np.array([3.8135955291021233, 39.889510090195238, 44.981835351176571, 2.0055189192561507])

        # Inter-individual variability matrix (OMEGA)
        self.OM = np.array([
            [0.10855133849022583, -0.0152445093837639736, -0.19698189309256298, 0.11914555131547180],
            [-0.0152445093837639736, 0.00301016276351715687, 0.023128525633882277, -0.006615975862019478],
            [-0.19698189309256298, 0.023128525633882277, 1.0420667710781930, -0.19223488058085114],
            [0.11914555131547180, -0.006615975862019478, -0.19223488058085114, 0.416]])

        # Residual variability matrix (SIGMA)
        self.SG = np.array([[0.14019731106615912 ** 2, 0.0],
                            [0.0, 1.8662549226759475 ** 2]])

    @staticmethod
    def Pred(TH, ETA, DATAi):

        V1 = TH[1] * np.exp(ETA[1])
        V2 = TH[2] * np.exp(ETA[2])
        Q = TH[3] * np.exp(ETA[3])
        K12 = Q / V1
        K21 = Q / V2

        Conc = pd.DataFrame([[0.0, 0.0]], columns=["C1", "C2"])
        IPRE = [0.0]
        inf = False

        pCLCR = min(DATAi.iloc[0]["CLCR"], 150)

        for i in range(1, len(DATAi)):
            # if i==2: break

            row = DATAi.iloc[i]
            prev = DATAi.iloc[i - 1]
            CLCR = min(row["CLCR"], 150) if not pd.isna(row["CLCR"]) else pCLCR
            pCLCR = CLCR

            TVCL = TH[0] * CLCR / 100
            CL = TVCL * np.exp(ETA[0])
            K10 = CL / V1

            AlpBe = K10 + K12 + K21
            AlmBe = K10 * K21
            Det4 = np.sqrt(AlpBe**2 - 4 * AlmBe)
            Alpha = (AlpBe + Det4) / 2
            Beta = (AlpBe - Det4) / 2
            Divisor = V1 * (Alpha - Beta)

            cTime = row["TIME"]
            dTime = cTime - prev["TIME"]

            if not pd.isna(prev["AMT"]) and prev["AMT"] > 0:
                pTime = prev["TIME"]
                pAMT = prev["AMT"]
                pRate = prev["RATE"]
                pDur = pAMT / pRate

                if dTime <= pDur:
                    infC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * dTime)) / Alpha
                    infC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * dTime)) / Beta
                    C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime) + infC1
                    C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime) + infC2
                    inf = True
                else:
                    eC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * pDur)) / Alpha
                    eC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * pDur)) / Beta
                    pC1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * pDur)
                    pC2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * pDur)
                    C1 = (eC1 + pC1) * np.exp(-Alpha * (dTime - pDur))
                    C2 = (eC2 + pC2) * np.exp(-Beta * (dTime - pDur))
                    inf = False
            elif inf:
                # pTime = prev["TIME"]
                # pAMT = prev["AMT"]
                # pRate = prev["RATE"]
                # pDur = pAMT / pRate
                if cTime <= pTime + pDur:
                    infC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * dTime)) / Alpha
                    infC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * dTime)) / Beta
                    C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime) + infC1
                    C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime) + infC2
                    inf = True
                else:
                    rDur = pTime + pDur - prev["TIME"]
                    eC1 = pRate * (Alpha - K21) / Divisor * (1 - np.exp(-Alpha * rDur)) / Alpha
                    eC2 = pRate * (K21 - Beta) / Divisor * (1 - np.exp(-Beta * rDur)) / Beta
                    pC1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * rDur)
                    pC2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * rDur)
                    C1 = (eC1 + pC1) * np.exp(-Alpha * (cTime - pTime - pDur))
                    C2 = (eC2 + pC2) * np.exp(-Beta * (cTime - pTime - pDur))
                    inf = False
            else:
                C1 = Conc.iloc[i - 1]["C1"] * np.exp(-Alpha * dTime)
                C2 = Conc.iloc[i - 1]["C2"] * np.exp(-Beta * dTime)
                inf = False

            Conc = pd.concat([Conc, pd.DataFrame([[C1, C2]], columns=["C1", "C2"])], ignore_index=True)
            IPRE.append(C1 + C2)

        return np.array(IPRE)