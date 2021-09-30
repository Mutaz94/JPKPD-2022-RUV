Wed Sep 29 11:25:40 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/B/dat67.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m67.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.35428831204        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1700E+02 -6.9342E+01  6.4095E+00 -1.0446E+02 -8.7724E+00  4.2235E+01 -1.4778E+01  4.2749E+00 -1.2105E+01  6.0903E+00
             6.4069E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.46698188458        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9388E-01  1.0775E+00  1.0439E+00  1.0678E+00  1.0575E+00  1.0458E+00  1.1021E+00  9.7640E-01  1.0208E+00  9.8497E-01
             9.7709E-01
 PARAMETER:  9.3866E-02  1.7466E-01  1.4300E-01  1.6563E-01  1.5593E-01  1.4482E-01  1.9723E-01  7.6112E-02  1.2061E-01  8.4860E-02
             7.6822E-02
 GRADIENT:   2.0305E-01 -3.1901E+00  2.6869E+00 -2.0221E+00  1.9333E+00  9.5585E+00 -7.1387E+00 -8.8912E-01  3.7837E+00 -6.1747E+00
            -6.6402E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.22877406956        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      353
 NPARAMETR:  9.9759E-01  1.0199E+00  1.1045E+00  1.1091E+00  1.0782E+00  1.0191E+00  1.2844E+00  9.1807E-01  9.3536E-01  1.0582E+00
             9.8822E-01
 PARAMETER:  9.7589E-02  1.1974E-01  1.9937E-01  2.0354E-01  1.7531E-01  1.1896E-01  3.5025E-01  1.4519E-02  3.3172E-02  1.5662E-01
             8.8145E-02
 GRADIENT:   9.7683E+00  1.4135E+00  1.6657E+00  4.1448E-01  7.7886E+00 -1.8715E-01 -2.3112E-01 -2.1351E+00 -1.7259E+00 -1.5794E+00
            -3.0639E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.52249022461        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.9096E-01  8.7266E-01  1.1382E+00  1.2049E+00  1.0127E+00  1.0171E+00  1.3974E+00  9.5088E-01  9.0920E-01  1.0244E+00
             9.9146E-01
 PARAMETER:  9.0916E-02 -3.6211E-02  2.2944E-01  2.8642E-01  1.1265E-01  1.1698E-01  4.3462E-01  4.9633E-02  4.8063E-03  1.2410E-01
             9.1418E-02
 GRADIENT:  -2.1170E+00  4.5172E+00  2.1512E+00  3.8626E+00 -5.1611E+00 -4.1010E-01  2.4474E-04 -1.5708E-01  8.4389E-01  6.6498E-01
             3.2659E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.59144977890        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.8867E-01  6.9808E-01  1.2833E+00  1.3215E+00  1.0120E+00  1.0142E+00  1.5782E+00  1.0622E+00  8.7005E-01  1.0478E+00
             9.9112E-01
 PARAMETER:  8.8604E-02 -2.5943E-01  3.4947E-01  3.7880E-01  1.1190E-01  1.1413E-01  5.5628E-01  1.6034E-01 -3.9203E-02  1.4667E-01
             9.1078E-02
 GRADIENT:  -1.7292E+00  5.0997E+00  3.0648E+00  5.8049E+00 -5.4546E+00 -5.7021E-01  3.7610E-01 -1.3954E-01  7.8345E-02  5.2469E-01
            -6.6345E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.63205526554        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.8641E-01  5.2615E-01  1.4579E+00  1.4402E+00  1.0260E+00  1.0117E+00  1.7880E+00  1.2149E+00  8.4147E-01  1.0808E+00
             9.9115E-01
 PARAMETER:  8.6319E-02 -5.4217E-01  4.7703E-01  4.6477E-01  1.2564E-01  1.1165E-01  6.8112E-01  2.9463E-01 -7.2605E-02  1.7774E-01
             9.1110E-02
 GRADIENT:  -6.5149E-01  6.7610E+00  3.7711E+00  1.1458E+01 -5.5315E+00 -5.3835E-01  9.4634E-01 -4.9354E-02 -6.0670E-01  4.5971E-01
            -5.0453E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.64891924511        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  9.8535E-01  4.3358E-01  1.5742E+00  1.5057E+00  1.0420E+00  1.0104E+00  1.9491E+00  1.3253E+00  8.2372E-01  1.1018E+00
             9.9095E-01
 PARAMETER:  8.5240E-02 -7.3569E-01  5.5376E-01  5.0923E-01  1.4115E-01  1.1037E-01  7.6739E-01  3.8167E-01 -9.3929E-02  1.9694E-01
             9.0913E-02
 GRADIENT:   3.3893E-01  7.2538E+00  3.7412E+00  1.6548E+01 -4.3766E+00 -5.0281E-01  9.9180E-01  1.4197E-01 -1.3058E+00  4.2520E-01
            -8.9466E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.67008862283        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.8468E-01  3.5604E-01  1.6673E+00  1.5596E+00  1.0531E+00  1.0097E+00  2.1430E+00  1.4101E+00  8.0752E-01  1.1142E+00
             9.9086E-01
 PARAMETER:  8.4560E-02 -9.3271E-01  6.1118E-01  5.4441E-01  1.5170E-01  1.0970E-01  8.6218E-01  4.4364E-01 -1.1379E-01  2.0813E-01
             9.0818E-02
 GRADIENT:   1.5477E+00  6.8780E+00  3.5406E+00  2.0244E+01 -2.5284E+00 -3.3119E-01  8.9494E-01 -5.8086E-02 -2.0753E+00  9.8381E-02
            -1.3039E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1664.72943540752        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  9.8392E-01  2.6875E-01  1.7642E+00  1.6186E+00  1.0616E+00  1.0092E+00  2.4595E+00  1.5008E+00  7.9028E-01  1.1222E+00
             9.9080E-01
 PARAMETER:  8.3790E-02 -1.2140E+00  6.6770E-01  5.8157E-01  1.5980E-01  1.0919E-01  9.9995E-01  5.0597E-01 -1.3536E-01  2.1528E-01
             9.0759E-02
 GRADIENT:   2.8338E+00  5.7002E+00  2.9329E+00  2.2200E+01 -2.0576E-01 -5.4861E-02  7.3848E-01 -2.9768E-01 -2.7256E+00 -3.7027E-01
            -1.6205E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1664.90132565994        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1585
 NPARAMETR:  9.8271E-01  1.6571E-01  1.8597E+00  1.6856E+00  1.0641E+00  1.0086E+00  3.1147E+00  1.5925E+00  7.7251E-01  1.1246E+00
             9.9101E-01
 PARAMETER:  8.2560E-02 -1.6975E+00  7.2042E-01  6.2211E-01  1.6209E-01  1.0852E-01  1.2361E+00  5.6530E-01 -1.5811E-01  2.1745E-01
             9.0972E-02
 GRADIENT:   3.6363E+00  3.6641E+00  1.9282E+00  2.0495E+01  1.9726E+00  2.3651E-01  4.8487E-01 -4.2715E-01 -2.9258E+00 -8.3635E-01
            -1.6607E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1665.22639069372        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1761
 NPARAMETR:  9.8065E-01  7.5147E-02  1.8707E+00  1.7400E+00  1.0410E+00  1.0074E+00  4.5728E+00  1.6104E+00  7.6101E-01  1.1118E+00
             9.9166E-01
 PARAMETER:  8.0463E-02 -2.4883E+00  7.2633E-01  6.5389E-01  1.4021E-01  1.0734E-01  1.6201E+00  5.7647E-01 -1.7311E-01  2.0601E-01
             9.1627E-02
 GRADIENT:   2.3014E+00  1.6466E+00  1.0271E+00  1.5879E+01  7.1621E-01  2.9021E-01  1.9991E-01 -1.8530E-01 -2.0861E+00 -7.3143E-01
            -9.6526E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1665.63244690306        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1943
 NPARAMETR:  9.7942E-01  2.6527E-02  1.8422E+00  1.7454E+00  1.0218E+00  1.0064E+00  6.5034E+00  1.5894E+00  7.5990E-01  1.1063E+00
             9.9241E-01
 PARAMETER:  7.9202E-02 -3.5296E+00  7.1094E-01  6.5698E-01  1.2154E-01  1.0638E-01  1.9723E+00  5.6333E-01 -1.7457E-01  2.0106E-01
             9.2386E-02
 GRADIENT:   1.7992E+00  8.5866E-02  3.8412E-02 -3.2760E+01  3.4148E+00  3.1118E-01 -2.2859E-02  2.7650E-01  4.9899E-01 -2.9520E-01
             2.7633E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.64509324359        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     2083
 NPARAMETR:  9.7857E-01  2.5458E-02  1.8404E+00  1.7462E+00  1.0173E+00  1.0056E+00  7.8767E+00  1.5858E+00  7.5866E-01  1.1075E+00
             9.9177E-01
 PARAMETER:  7.8332E-02 -3.5707E+00  7.0996E-01  6.5743E-01  1.1719E-01  1.0559E-01  2.1639E+00  5.6110E-01 -1.7621E-01  2.0211E-01
             9.1733E-02
 GRADIENT:  -4.7090E-03  2.2660E-01  1.5612E+00 -3.1964E+01 -1.2173E-02 -7.5079E-05  1.6602E-01  3.2225E-01  3.2980E-01  3.7650E-01
             8.4539E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1665.64797808638        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2258
 NPARAMETR:  9.7501E-01  2.0467E-02  1.8243E+00  1.7501E+00  1.0134E+00  1.0031E+00  8.8197E+00  1.5736E+00  7.5668E-01  1.1012E+00
             9.9176E-01
 PARAMETER:  7.4694E-02 -3.7889E+00  7.0122E-01  6.5967E-01  1.1336E-01  1.0311E-01  2.2770E+00  5.5337E-01 -1.7881E-01  1.9642E-01
             9.1730E-02
 GRADIENT:  -7.7972E+00  2.6662E-01  9.9381E-02 -2.9086E+01  2.2030E+00 -9.8727E-01  2.5158E-01  2.4560E-01 -4.7491E-01 -2.6989E-01
             9.6939E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1665.68279560883        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2445             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7915E-01  1.7909E-02  1.8177E+00  1.7505E+00  1.0078E+00  1.0059E+00  8.5299E+00  1.5654E+00  7.5784E-01  1.1025E+00
             9.9162E-01
 PARAMETER:  7.8930E-02 -3.9224E+00  6.9759E-01  6.5992E-01  1.0778E-01  1.0593E-01  2.2436E+00  5.4811E-01 -1.7728E-01  1.9762E-01
             9.1589E-02
 GRADIENT:   4.0570E+02  2.8099E+00  1.1190E+01  1.2284E+03  4.8921E+00  5.2548E+01  1.0903E+00  2.3424E+00  2.8436E+01  2.2216E+00
             9.0506E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1665.69409273817        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2629
 NPARAMETR:  9.7891E-01  1.1436E-02  1.8140E+00  1.7557E+00  1.0077E+00  1.0058E+00  9.4258E+00  1.5625E+00  7.5758E-01  1.1000E+00
             9.9153E-01
 PARAMETER:  7.8684E-02 -4.3710E+00  6.9554E-01  6.6289E-01  1.0769E-01  1.0577E-01  2.3434E+00  5.4631E-01 -1.7762E-01  1.9528E-01
             9.1494E-02
 GRADIENT:   1.3563E+00  7.4160E-02  1.8542E-02 -2.8014E+01  1.7426E+00  1.7398E-01 -1.1957E-02  1.3286E-01  3.8098E-01 -5.1744E-02
             4.4917E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1665.70082167648        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2796             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7853E-01  1.0000E-02  1.8136E+00  1.7560E+00  1.0063E+00  1.0055E+00  1.1738E+01  1.5615E+00  7.5709E-01  1.1001E+00
             9.9153E-01
 PARAMETER:  7.8300E-02 -4.5492E+00  6.9529E-01  6.6306E-01  1.0628E-01  1.0547E-01  2.5628E+00  5.4562E-01 -1.7827E-01  1.9536E-01
             9.1490E-02
 GRADIENT:   4.0460E+02  2.1374E-03  1.0273E+01  1.2429E+03  6.5023E+00  5.2216E+01  7.4290E-01  2.1956E+00  2.9008E+01  1.8361E+00
             8.3340E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1665.70710942087        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     2958             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7895E-01  1.0000E-02  1.8046E+00  1.7558E+00  1.0025E+00  1.0058E+00  1.2044E+01  1.5541E+00  7.5678E-01  1.0972E+00
             9.9148E-01
 PARAMETER:  7.8720E-02 -4.5531E+00  6.9032E-01  6.6290E-01  1.0249E-01  1.0576E-01  2.5886E+00  5.4089E-01 -1.7869E-01  1.9275E-01
             9.1446E-02
 GRADIENT:   4.0541E+02  0.0000E+00  1.0529E+01  1.2425E+03  5.4795E+00  5.2399E+01  8.0105E-01  2.2049E+00  2.8907E+01  1.8319E+00
             8.5517E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1665.71002274759        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3147
 NPARAMETR:  9.7895E-01  1.0000E-02  1.7989E+00  1.7555E+00  1.0007E+00  1.0058E+00  1.2005E+01  1.5488E+00  7.5672E-01  1.0964E+00
             9.9131E-01
 PARAMETER:  7.8729E-02 -4.5531E+00  6.8720E-01  6.6273E-01  1.0069E-01  1.0575E-01  2.5853E+00  5.3746E-01 -1.7876E-01  1.9199E-01
             9.1270E-02
 GRADIENT:   1.4853E+00  0.0000E+00  6.3268E-01 -2.8562E+01 -3.5332E-02  1.7788E-01  1.6954E-02  1.5453E-01  1.1197E-01  1.6567E-01
             1.8141E-03

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1665.71288337998        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3339             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7894E-01  1.0000E-02  1.7920E+00  1.7551E+00  9.9926E-01  1.0058E+00  1.2020E+01  1.5432E+00  7.5687E-01  1.0935E+00
             9.9128E-01
 PARAMETER:  7.8716E-02 -4.5531E+00  6.8335E-01  6.6254E-01  9.9259E-02  1.0574E-01  2.5865E+00  5.3384E-01 -1.7856E-01  1.8942E-01
             9.1239E-02
 GRADIENT:   4.0537E+02  0.0000E+00  1.0048E+01  1.2414E+03  6.1084E+00  5.2365E+01  7.9659E-01  2.1068E+00  2.8843E+01  1.5520E+00
             7.9181E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1665.71410404366        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3527
 NPARAMETR:  9.7894E-01  1.0000E-02  1.7863E+00  1.7548E+00  9.9831E-01  1.0057E+00  1.2018E+01  1.5385E+00  7.5695E-01  1.0924E+00
             9.9118E-01
 PARAMETER:  7.8710E-02 -4.5531E+00  6.8013E-01  6.6238E-01  9.8305E-02  1.0573E-01  2.5864E+00  5.3082E-01 -1.7845E-01  1.8839E-01
             9.1145E-02
 GRADIENT:   1.6508E+00  0.0000E+00 -3.1433E-01 -2.8833E+01  1.5798E+00  1.8110E-01  1.7024E-02  1.1203E-01  1.2342E-01 -2.2865E-01
            -3.3237E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1665.71657714871        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3719             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7893E-01  1.0000E-02  1.7832E+00  1.7546E+00  9.9571E-01  1.0057E+00  1.2027E+01  1.5342E+00  7.5707E-01  1.0931E+00
             9.9121E-01
 PARAMETER:  7.8701E-02 -4.5531E+00  6.7839E-01  6.6225E-01  9.5700E-02  1.0573E-01  2.5871E+00  5.2801E-01 -1.7830E-01  1.8902E-01
             9.1167E-02
 GRADIENT:   4.0535E+02  0.0000E+00  1.0335E+01  1.2404E+03  5.2395E+00  5.2348E+01  7.9776E-01  2.0715E+00  2.8836E+01  1.7814E+00
             8.1354E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1665.71773865968        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3907
 NPARAMETR:  9.7892E-01  1.0000E-02  1.7794E+00  1.7544E+00  9.9444E-01  1.0057E+00  1.2028E+01  1.5306E+00  7.5714E-01  1.0925E+00
             9.9117E-01
 PARAMETER:  7.8695E-02 -4.5531E+00  6.7626E-01  6.6214E-01  9.4429E-02  1.0572E-01  2.5873E+00  5.2567E-01 -1.7821E-01  1.8848E-01
             9.1127E-02
 GRADIENT:   1.5277E+00  0.0000E+00  4.2979E-01 -2.8650E+01 -1.9145E-01  1.7377E-01  1.7418E-02  9.8550E-02  1.4222E-01  1.4553E-01
             1.2274E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1665.71857021466        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     4087
 NPARAMETR:  9.7891E-01  1.0000E-02  1.7717E+00  1.7541E+00  9.9384E-01  1.0057E+00  1.2026E+01  1.5251E+00  7.5723E-01  1.0896E+00
             9.9107E-01
 PARAMETER:  7.8691E-02 -4.5531E+00  6.7465E-01  6.6205E-01  9.3610E-02  1.0572E-01  2.5874E+00  5.2402E-01 -1.7814E-01  1.8772E-01
             9.1108E-02
 GRADIENT:   1.2248E-03  0.0000E+00  3.8952E-01  1.6478E-01 -7.6062E-02  3.3943E-04  5.0942E-05  7.4434E-02 -7.4352E-03  9.3598E-02
             1.2024E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4087
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0904E-04  3.9892E-04 -3.6457E-02 -7.3466E-03 -5.0194E-02
 SE:             2.9843E-02  1.8059E-03  1.9169E-02  2.9237E-02  1.9578E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9441E-01  8.2517E-01  5.7187E-02  8.0160E-01  1.0355E-02

 ETASHRINKSD(%)  2.1664E-02  9.3950E+01  3.5782E+01  2.0536E+00  3.4411E+01
 ETASHRINKVR(%)  4.3324E-02  9.9634E+01  5.8760E+01  4.0650E+00  5.6980E+01
 EBVSHRINKSD(%)  4.0500E-01  9.4149E+01  3.9050E+01  2.3645E+00  3.0966E+01
 EBVSHRINKVR(%)  8.0836E-01  9.9658E+01  6.2850E+01  4.6732E+00  5.2343E+01
 RELATIVEINF(%)  9.3566E+01  9.1134E-03  1.0575E+01  2.9715E+00  7.2521E+00
 EPSSHRINKSD(%)  4.5358E+01
 EPSSHRINKVR(%)  7.0143E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.7185702146589     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.56774365092076     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    61.42
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     6.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.719       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  1.78E+00  1.75E+00  9.94E-01  1.01E+00  1.20E+01  1.53E+00  7.57E-01  1.09E+00  9.91E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.96E-02  0.00E+00  2.94E-01  4.62E-02  9.46E-02  6.49E-02  3.80E-03  3.29E-01  6.16E-02  1.99E-01  7.06E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.75E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -7.71E-05  0.00E+00  8.62E-02
 
 TH 4
+       -1.53E-05  0.00E+00  6.63E-03  2.13E-03
 
 TH 5
+       -2.89E-05  0.00E+00  2.40E-02  2.20E-03  8.94E-03
 
 TH 6
+        7.35E-05  0.00E+00  2.08E-03 -4.70E-05  6.58E-04  4.22E-03
 
 TH 7
+        3.06E-06  0.00E+00  2.86E-04  1.59E-05  9.26E-05 -1.28E-05  1.45E-05
 
 TH 8
+       -6.12E-05  0.00E+00  6.86E-02  4.11E-03  1.88E-02  3.07E-03  3.46E-04  1.08E-01
 
 TH 9
+       -1.67E-06  0.00E+00 -1.72E-03 -8.21E-04 -8.25E-04 -2.25E-04  1.04E-04 -4.11E-04  3.79E-03
 
 TH10
+       -3.24E-04  0.00E+00  8.14E-03  9.62E-04  5.18E-03 -1.46E-03 -7.37E-05 -1.02E-02 -1.46E-03  3.95E-02
 
 TH11
+       -2.23E-04  0.00E+00  2.73E-03  4.60E-04  7.12E-04 -2.20E-04  7.66E-06 -1.56E-03 -3.52E-04  1.51E-03  4.98E-03
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.96E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.88E-03  0.00E+00  2.94E-01
 
 TH 4
+       -1.12E-02  0.00E+00  4.89E-01  4.62E-02
 
 TH 5
+       -1.03E-02  0.00E+00  8.65E-01  5.03E-01  9.46E-02
 
 TH 6
+        3.83E-02  0.00E+00  1.09E-01 -1.57E-02  1.07E-01  6.49E-02
 
 TH 7
+        2.72E-02  0.00E+00  2.56E-01  9.05E-02  2.57E-01 -5.19E-02  3.80E-03
 
 TH 8
+       -6.28E-03  0.00E+00  7.10E-01  2.71E-01  6.04E-01  1.44E-01  2.77E-01  3.29E-01
 
 TH 9
+       -9.19E-04  0.00E+00 -9.53E-02 -2.89E-01 -1.42E-01 -5.63E-02  4.45E-01 -2.03E-02  6.16E-02
 
 TH10
+       -5.52E-02  0.00E+00  1.39E-01  1.05E-01  2.76E-01 -1.13E-01 -9.74E-02 -1.55E-01 -1.19E-01  1.99E-01
 
 TH11
+       -1.07E-01  0.00E+00  1.32E-01  1.41E-01  1.07E-01 -4.81E-02  2.85E-02 -6.72E-02 -8.11E-02  1.08E-01  7.06E-02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.16E+03
 
 TH 2
+        1.70E-14 -2.05E-29
 
 TH 3
+       -4.02E+00  4.80E-14  3.62E+01
 
 TH 4
+        5.42E+00  1.05E-14  2.19E+00  6.46E+02
 
 TH 5
+       -1.26E+01 -1.14E-13 -1.26E+02 -8.65E+01  4.66E+02
 
 TH 6
+       -1.46E+01  7.67E-14  1.30E+01  2.20E+01 -3.67E+01  2.34E+02
 
 TH 7
+        1.99E-03  4.54E-18  1.09E-03  1.23E-02 -2.28E-03  1.66E-03  5.78E-06
 
 TH 8
+        4.46E+00 -1.21E-14  1.17E-01 -6.06E+00 -4.22E+00 -2.37E+00 -1.72E-04  1.11E+00
 
 TH 9
+        5.51E+00  2.63E-14  2.12E+00  1.08E+02  8.87E-01  9.72E+00  4.04E-02 -1.59E+00  2.83E+02
 
 TH10
+        7.74E+00  3.37E-15  1.72E+01 -5.92E-01 -6.69E+01 -3.19E+00 -5.60E-04  1.65E+00 -6.76E+00  1.09E+01
 
 TH11
+        5.58E+01 -1.83E-13 -2.32E+01 -1.97E+01  2.48E+01  1.36E+01  2.41E-03  1.41E+01  1.72E+01  7.80E+00  2.16E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.55E+00  0.00E+00  5.84E+01
 
 TH 4
+       -2.52E+01  0.00E+00 -1.88E+00  6.01E+02
 
 TH 5
+        2.77E+01  0.00E+00 -1.25E+02 -1.44E+00  4.82E+02
 
 TH 6
+        1.50E+01  0.00E+00 -2.10E+00 -2.24E+01  2.29E+01  1.58E+02
 
 TH 7
+        3.29E-03  0.00E+00  2.24E-03 -2.24E-02 -6.74E-04 -2.68E-03  1.98E-05
 
 TH 8
+       -5.84E+00  0.00E+00 -1.74E+01 -2.25E+01 -9.85E+00  6.85E-01  2.64E-03  2.04E+01
 
 TH 9
+        4.14E+00  0.00E+00  1.79E+01 -1.28E+02 -1.89E+01 -1.20E+01  7.66E-02  4.23E+00  3.78E+02
 
 TH10
+       -1.58E+01  0.00E+00  3.49E+00 -1.38E+01 -1.08E+02 -1.62E+01 -4.14E-03  1.33E+01 -1.38E+01  9.26E+01
 
 TH11
+       -7.33E+01  0.00E+00  1.81E+00  5.14E+00 -4.47E+01 -7.84E+00  1.72E-03  7.88E+00  5.16E+00  2.21E+01  2.06E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       68.461
Stop Time:
Wed Sep 29 11:26:50 CDT 2021
