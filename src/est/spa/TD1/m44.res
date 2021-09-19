Sat Sep 18 14:04:04 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat44.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m44.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1579.70768563624        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9368E+02  3.2330E+01  1.5329E+01  4.6000E+01  1.7464E+01 -2.0642E+01  1.5129E+00 -1.8556E+01  1.1096E+01 -9.3345E+00
            -3.4507E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1584.91461511636        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9523E-01  9.7232E-01  8.9499E-01  1.0098E+00  9.0774E-01  1.0417E+00  9.8386E-01  1.1976E+00  8.8860E-01  9.7126E-01
             1.0657E+00
 PARAMETER:  9.5214E-02  7.1927E-02 -1.0940E-02  1.0972E-01  3.1973E-03  1.4082E-01  8.3731E-02  2.8032E-01 -1.8104E-02  7.0837E-02
             1.6361E-01
 GRADIENT:   1.6361E+02  3.3675E+01  2.3820E+00  4.4526E+01 -2.1606E+01 -3.1829E-01 -4.3184E+00  1.3088E+00 -6.7530E+00  9.3028E+00
            -2.4078E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1587.80540226940        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.6942E-01  9.3089E-01  8.9323E-01  1.0239E+00  8.9737E-01  1.0343E+00  1.2938E+00  1.2470E+00  8.0043E-01  8.3672E-01
             1.0684E+00
 PARAMETER:  6.8947E-02  2.8386E-02 -1.2911E-02  1.2359E-01 -8.2918E-03  1.3371E-01  3.5758E-01  3.2070E-01 -1.2260E-01 -7.8263E-02
             1.6614E-01
 GRADIENT:   1.1192E+02  2.2818E+01 -2.2826E+00  2.7856E+01  4.3124E+00  2.0448E+00  9.9374E+00  2.5545E+00 -1.9309E+00 -4.3409E+00
            -2.9669E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1590.27490170009        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.2782E-01  8.8262E-01  9.6893E-01  1.0462E+00  9.1593E-01  1.0152E+00  1.1872E+00  1.2372E+00  8.3601E-01  9.0430E-01
             1.0669E+00
 PARAMETER:  2.5086E-02 -2.4855E-02  6.8435E-02  1.4514E-01  1.2189E-02  1.1511E-01  2.7163E-01  3.1285E-01 -7.9115E-02 -5.9621E-04
             1.6474E-01
 GRADIENT:   2.3285E+01  7.9367E+00  1.4118E+00  1.0202E+01 -2.0945E+00 -2.3849E+00  1.3778E+00  3.7569E-01 -1.7114E+00 -5.3765E-01
            -1.8315E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1591.22700840936        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      389
 NPARAMETR:  9.2851E-01  6.7313E-01  1.2303E+00  1.1900E+00  9.4911E-01  1.0324E+00  1.2196E+00  1.4130E+00  8.2307E-01  9.6887E-01
             1.0713E+00
 PARAMETER:  2.5828E-02 -2.9581E-01  3.0730E-01  2.7393E-01  4.7772E-02  1.3188E-01  2.9849E-01  4.4569E-01 -9.4715E-02  6.8371E-02
             1.6891E-01
 GRADIENT:  -3.4906E+00  4.4931E+00 -8.6559E-01  4.2910E+00 -1.4309E-01  6.5580E-02  8.3353E-01  1.3307E+00 -1.0183E+00 -1.1040E+00
            -3.3839E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1592.68938933634        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      566
 NPARAMETR:  9.2820E-01  4.1151E-01  1.5355E+00  1.3651E+00  9.7069E-01  1.0265E+00  9.1178E-01  1.5832E+00  7.9839E-01  1.0601E+00
             1.0731E+00
 PARAMETER:  2.5487E-02 -7.8793E-01  5.2887E-01  4.1125E-01  7.0252E-02  1.2615E-01  7.6404E-03  5.5947E-01 -1.2516E-01  1.5835E-01
             1.7051E-01
 GRADIENT:   4.0894E+00  4.5626E+00  3.0054E+00  1.6333E+01 -6.1487E+00 -7.3014E-01  9.4634E-01 -3.6966E-01  1.9184E+00  1.5992E+00
             3.4173E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1593.06231549527        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      742
 NPARAMETR:  9.2487E-01  3.0628E-01  1.5534E+00  1.4221E+00  9.4949E-01  1.0267E+00  6.1457E-01  1.6067E+00  7.6942E-01  1.0344E+00
             1.0715E+00
 PARAMETER:  2.1896E-02 -1.0833E+00  5.4044E-01  4.5212E-01  4.8171E-02  1.2631E-01 -3.8682E-01  5.7417E-01 -1.6212E-01  1.3387E-01
             1.6904E-01
 GRADIENT:  -5.4020E-02 -5.0889E-01 -2.5725E-01 -9.7372E-01  4.1075E-01 -1.9659E-01  2.6862E-01  6.9760E-02  2.2366E-01 -1.0096E-01
             2.5223E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1593.28371389541        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      922
 NPARAMETR:  9.2308E-01  4.2802E-01  1.4864E+00  1.3515E+00  9.6626E-01  1.0257E+00  1.3228E-01  1.5596E+00  8.2022E-01  1.0549E+00
             1.0700E+00
 PARAMETER:  1.9956E-02 -7.4859E-01  4.9637E-01  4.0119E-01  6.5677E-02  1.2533E-01 -1.9229E+00  5.4444E-01 -9.8184E-02  1.5343E-01
             1.6770E-01
 GRADIENT:  -8.1865E+00  3.6357E+00 -1.2357E+00  1.8919E+01 -1.6990E+00 -1.2551E+00  3.4958E-02  5.6807E-01  8.8863E-02  8.1500E-01
            -4.5916E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1593.40773693864        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1097
 NPARAMETR:  9.2735E-01  4.9221E-01  1.4812E+00  1.3027E+00  9.8506E-01  1.0295E+00  2.0813E-02  1.5635E+00  8.4628E-01  1.0614E+00
             1.0719E+00
 PARAMETER:  2.4579E-02 -6.0885E-01  4.9282E-01  3.6445E-01  8.4947E-02  1.2907E-01 -3.7722E+00  5.4690E-01 -6.6904E-02  1.5957E-01
             1.6941E-01
 GRADIENT:  -1.5720E-01 -5.8138E-02 -1.2235E-01 -3.3443E-02  2.7322E-01  8.6701E-03  1.3331E-03 -4.9311E-03 -1.0629E-01 -8.4559E-02
            -1.9190E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1593.40870211274        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  9.2757E-01  5.0496E-01  1.4780E+00  1.2945E+00  9.8797E-01  1.0296E+00  1.4966E-02  1.5644E+00  8.5213E-01  1.0638E+00
             1.0719E+00
 PARAMETER:  2.4810E-02 -5.8328E-01  4.9066E-01  3.5811E-01  8.7894E-02  1.2916E-01 -4.1020E+00  5.4751E-01 -6.0018E-02  1.6181E-01
             1.6947E-01
 GRADIENT:  -1.7254E-02 -3.5458E-02 -4.9274E-03 -1.4857E-01  9.6163E-03 -1.2919E-02  7.3904E-04  4.7433E-03  3.3701E-02  3.0405E-03
             7.4298E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1593.40870602762        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1448
 NPARAMETR:  9.2757E-01  5.0381E-01  1.4786E+00  1.2953E+00  9.8780E-01  1.0296E+00  1.5062E-02  1.5645E+00  8.5153E-01  1.0637E+00
             1.0719E+00
 PARAMETER:  2.4811E-02 -5.8556E-01  4.9108E-01  3.5874E-01  8.7727E-02  1.2921E-01 -4.0956E+00  5.4759E-01 -6.0718E-02  1.6171E-01
             1.6945E-01
 GRADIENT:   1.4200E-02  5.4491E-03 -4.0459E-03 -1.9395E-03  1.8680E-02  1.1608E-02  7.4327E-04 -4.4337E-04 -3.9825E-03 -3.1865E-03
            -3.4452E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1593.40882572501        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1628
 NPARAMETR:  9.2756E-01  5.0059E-01  1.4801E+00  1.2975E+00  9.8725E-01  1.0297E+00  1.0000E-02  1.5648E+00  8.4994E-01  1.0633E+00
             1.0719E+00
 PARAMETER:  2.4805E-02 -5.9197E-01  4.9210E-01  3.6046E-01  8.7171E-02  1.2930E-01 -4.5340E+00  5.4773E-01 -6.2586E-02  1.6140E-01
             1.6943E-01
 GRADIENT:   8.8068E-02  8.5148E-02 -2.5477E-04  2.9931E-01  2.8695E-02  5.9938E-02  0.0000E+00 -1.0737E-02 -8.0293E-02 -1.4656E-02
            -2.4452E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1593.40890735971        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1794
 NPARAMETR:  9.2756E-01  5.0293E-01  1.4787E+00  1.2959E+00  9.8752E-01  1.0296E+00  1.0000E-02  1.5643E+00  8.5113E-01  1.0635E+00
             1.0719E+00
 PARAMETER:  2.4798E-02 -5.8730E-01  4.9114E-01  3.5919E-01  8.7441E-02  1.2918E-01 -4.5221E+00  5.4744E-01 -6.1195E-02  1.6155E-01
             1.6944E-01
 GRADIENT:   1.1136E-02  1.9734E-02  3.5103E-03  4.1890E-02 -9.5113E-03  3.3875E-03  0.0000E+00 -3.1169E-03 -1.7746E-02 -1.5033E-03
            -6.7227E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1794
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8577E-04 -2.5215E-04 -4.0879E-02 -5.2013E-03 -4.4091E-02
 SE:             2.9845E-02  9.8128E-05  1.8556E-02  2.9086E-02  2.0496E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7899E-01  1.0181E-02  2.7598E-02  8.5808E-01  3.1459E-02

 ETASHRINKSD(%)  1.6849E-02  9.9671E+01  3.7834E+01  2.5568E+00  3.1336E+01
 ETASHRINKVR(%)  3.3695E-02  9.9999E+01  6.1354E+01  5.0483E+00  5.2853E+01
 EBVSHRINKSD(%)  4.5802E-01  9.9683E+01  4.2281E+01  3.0418E+00  2.7576E+01
 EBVSHRINKVR(%)  9.1394E-01  9.9999E+01  6.6685E+01  5.9910E+00  4.7547E+01
 RELATIVEINF(%)  9.8246E+01  6.2123E-05  9.6151E+00  6.7422E+00  1.2263E+01
 EPSSHRINKSD(%)  4.5219E+01
 EPSSHRINKVR(%)  6.9990E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1593.4089073597072     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -858.25808079596902     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.15
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1593.409       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.28E-01  5.03E-01  1.48E+00  1.30E+00  9.88E-01  1.03E+00  1.00E-02  1.56E+00  8.51E-01  1.06E+00  1.07E+00
 


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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.20E+03
 
 TH 2
+       -2.37E+01  4.68E+02
 
 TH 3
+        1.81E+00  4.77E+01  7.46E+01
 
 TH 4
+       -1.27E+01  5.60E+02 -1.27E+01  8.59E+02
 
 TH 5
+        2.83E+00 -1.92E+02 -1.48E+02 -5.17E+01  5.45E+02
 
 TH 6
+       -9.61E-01 -3.60E+00  5.69E-01 -3.41E+00 -1.95E+00  1.83E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        6.26E-02 -1.01E+01 -2.14E+01 -3.12E+00 -4.31E+00  5.12E-02  0.00E+00  2.25E+01
 
 TH 9
+       -3.83E+00 -1.04E+02  4.85E+00  1.66E+00 -3.25E-01  3.60E-03  0.00E+00  5.99E-02  2.38E+02
 
 TH10
+        1.94E-01  7.28E+00 -2.58E+00  4.05E-01 -7.06E+01  7.96E-01  0.00E+00  1.07E+01 -7.27E-01  5.83E+01
 
 TH11
+       -9.94E+00 -1.51E+01 -7.06E+00 -1.02E+01 -4.46E+00  3.49E+00  0.00E+00  4.61E+00  1.01E+01  1.07E+01  1.77E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       26.986
Stop Time:
Sat Sep 18 14:04:32 CDT 2021
