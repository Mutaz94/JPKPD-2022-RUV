Sat Sep 25 13:40:20 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat63.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.09141999874        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.8675E+01 -1.0374E+02 -5.4669E+01 -1.0700E+02  3.3867E+01 -2.2475E+01 -9.6936E+00  1.6405E+01 -1.0405E+01  3.1200E+01
             2.7737E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1725.51434206394        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  1.0453E+00  1.0827E+00  1.1902E+00  1.0383E+00  1.0676E+00  1.0412E+00  1.0412E+00  9.1230E-01  1.0378E+00  8.5465E-01
             9.5082E-01
 PARAMETER:  1.4426E-01  1.7949E-01  2.7408E-01  1.3756E-01  1.6538E-01  1.4039E-01  1.4037E-01  8.2124E-03  1.3711E-01 -5.7065E-02
             4.9566E-02
 GRADIENT:   6.6471E+01  7.2532E+00  8.8856E-01  7.3150E+00  3.6971E+00  5.1193E-01 -2.1219E+00  1.8297E+00 -3.5228E-01 -3.0864E+00
            -1.6381E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1726.08957490320        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0391E+00  9.6554E-01  1.0859E+00  1.1128E+00  9.8058E-01  1.0462E+00  1.2993E+00  6.4044E-01  9.3754E-01  8.2836E-01
             9.4022E-01
 PARAMETER:  1.3837E-01  6.4933E-02  1.8238E-01  2.0685E-01  8.0393E-02  1.4518E-01  3.6185E-01 -3.4559E-01  3.5502E-02 -8.8306E-02
             3.8362E-02
 GRADIENT:   5.2754E+01  1.0099E+01 -1.9574E+00  2.5877E+01  7.9131E+00  2.5038E+00  5.7988E+00  3.2635E-01 -1.8370E+00  2.0040E+00
            -4.6916E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1726.56354468243        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      317
 NPARAMETR:  1.0479E+00  9.2682E-01  9.8297E-01  1.1232E+00  9.1178E-01  1.0595E+00  1.3154E+00  5.1354E-01  9.2750E-01  7.5349E-01
             9.5239E-01
 PARAMETER:  1.4683E-01  2.4000E-02  8.2827E-02  2.1614E-01  7.6454E-03  1.5777E-01  3.7411E-01 -5.6642E-01  2.4738E-02 -1.8304E-01
             5.1215E-02
 GRADIENT:   1.5393E+00  1.1501E+00 -7.7912E-01  1.1410E+00  3.5764E-01  4.1888E-01  2.2197E-02  3.8875E-01 -8.8452E-02  2.3208E-02
             3.8932E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.61384524522        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0445E+00  8.1080E-01  1.0351E+00  1.1989E+00  8.8973E-01  1.0570E+00  1.4481E+00  5.4615E-01  8.9186E-01  7.5605E-01
             9.5422E-01
 PARAMETER:  1.4354E-01 -1.0974E-01  1.3452E-01  2.8141E-01 -1.6836E-02  1.5542E-01  4.7023E-01 -5.0486E-01 -1.4450E-02 -1.7964E-01
             5.3138E-02
 GRADIENT:  -2.5787E+00  3.4071E+00  1.0015E+00  4.4298E+00 -2.6830E+00 -7.6420E-02  7.0993E-01  5.3803E-02  4.7013E-03 -8.0706E-02
             1.0584E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1726.73482573675        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.0448E+00  6.4575E-01  1.2933E+00  1.3169E+00  9.4122E-01  1.0535E+00  1.5213E+00  8.2404E-01  8.7175E-01  8.3979E-01
             9.4729E-01
 PARAMETER:  1.4381E-01 -3.3734E-01  3.5722E-01  3.7525E-01  3.9426E-02  1.5213E-01  5.1958E-01 -9.3533E-02 -3.7251E-02 -7.4601E-02
             4.5854E-02
 GRADIENT:   4.2254E+00  2.8525E+00  1.7888E+00  3.2620E+00 -4.9163E+00  3.3698E-02 -2.6720E-01  5.3198E-01  1.3589E-01  9.9901E-01
            -1.3841E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1727.01494695209        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      853
 NPARAMETR:  1.0368E+00  4.3583E-01  1.4255E+00  1.4551E+00  9.2545E-01  1.0503E+00  1.9069E+00  9.2055E-01  8.1727E-01  8.4322E-01
             9.5335E-01
 PARAMETER:  1.3615E-01 -7.3050E-01  4.5454E-01  4.7505E-01  2.2520E-02  1.4904E-01  7.4547E-01  1.7218E-02 -1.0179E-01 -7.0526E-02
             5.2222E-02
 GRADIENT:  -4.9636E+00  3.6317E+00  2.7413E+00  1.1815E+01 -3.4548E+00  1.0836E-01  2.9949E-01 -5.6216E-01 -8.5533E-01 -7.7765E-01
             8.4054E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1727.24427331336        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1029
 NPARAMETR:  1.0375E+00  2.7682E-01  1.5103E+00  1.5539E+00  9.1056E-01  1.0472E+00  2.3448E+00  1.0103E+00  7.8849E-01  8.6554E-01
             9.4851E-01
 PARAMETER:  1.3680E-01 -1.1844E+00  5.1230E-01  5.4079E-01  6.3054E-03  1.4609E-01  9.5222E-01  1.1026E-01 -1.3764E-01 -4.4401E-02
             4.7137E-02
 GRADIENT:   1.4327E+00  1.8044E+00  6.1083E-01  9.8921E+00 -4.9012E+00 -3.8062E-02  1.6577E-01  3.6308E-01 -1.8439E-01  1.2194E+00
            -2.5606E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1727.47653793873        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  1.0349E+00  1.3123E-01  1.7546E+00  1.6528E+00  9.4501E-01  1.0442E+00  3.0277E+00  1.2267E+00  7.6092E-01  8.8307E-01
             9.4883E-01
 PARAMETER:  1.3433E-01 -1.9308E+00  6.6222E-01  6.0247E-01  4.3441E-02  1.4323E-01  1.2078E+00  3.0433E-01 -1.7323E-01 -2.4356E-02
             4.7476E-02
 GRADIENT:   9.2125E-01  5.5144E-01  1.5923E+00  6.1469E+00 -3.0095E+00 -2.4662E-01 -4.0281E-01  1.1392E-01 -7.2150E-01 -1.8930E-01
            -5.2004E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1727.66855336214        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1381
 NPARAMETR:  1.0326E+00  4.6956E-02  1.7791E+00  1.7056E+00  9.3271E-01  1.0435E+00  5.3166E+00  1.2500E+00  7.4600E-01  8.8663E-01
             9.4945E-01
 PARAMETER:  1.3207E-01 -2.9585E+00  6.7609E-01  6.3393E-01  3.0336E-02  1.4255E-01  1.7708E+00  3.2318E-01 -1.9304E-01 -2.0327E-02
             4.8125E-02
 GRADIENT:  -9.1811E-01  2.3994E-01 -5.1908E-01  8.2167E+00 -1.5361E+00 -3.7404E-02 -2.6144E-02  4.2111E-01 -8.4655E-02  6.3789E-01
             2.3166E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1727.73050791367        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1556
 NPARAMETR:  1.0327E+00  1.6452E-02  1.8896E+00  1.7278E+00  9.5210E-01  1.0432E+00  9.0914E+00  1.3311E+00  7.3924E-01  8.8955E-01
             9.4960E-01
 PARAMETER:  1.3215E-01 -4.0073E+00  7.3635E-01  6.4685E-01  5.0919E-02  1.4229E-01  2.3073E+00  3.8603E-01 -2.0214E-01 -1.7041E-02
             4.8282E-02
 GRADIENT:  -1.4117E-01  9.8253E-02  1.6802E+00  4.3979E+00 -2.0788E+00  5.7975E-03  3.6956E-02 -2.1912E-01 -1.0725E-01 -2.3120E-01
            -1.4427E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1727.75629331295        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1733
 NPARAMETR:  1.0323E+00  1.0037E-02  1.8603E+00  1.7293E+00  9.4471E-01  1.0431E+00  1.1746E+01  1.3127E+00  7.3847E-01  8.8900E-01
             9.4969E-01
 PARAMETER:  1.3178E-01 -4.5015E+00  7.2074E-01  6.4771E-01  4.3123E-02  1.4221E-01  2.5635E+00  3.7208E-01 -2.0317E-01 -1.7658E-02
             4.8378E-02
 GRADIENT:  -5.4256E-01  2.0626E-01  1.9225E-01  1.8332E+00 -4.8522E-01  3.0911E-02  3.1185E-01  9.0506E-02 -2.6582E-01 -2.4573E-02
             1.4424E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1727.75807717301        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1918
 NPARAMETR:  1.0324E+00  1.0012E-02  1.8568E+00  1.7281E+00  9.4445E-01  1.0430E+00  1.1729E+01  1.3090E+00  7.3887E-01  8.8885E-01
             9.4902E-01
 PARAMETER:  1.3191E-01 -4.5040E+00  7.1886E-01  6.4702E-01  4.2848E-02  1.4209E-01  2.5621E+00  3.6923E-01 -2.0264E-01 -1.7822E-02
             4.7671E-02
 GRADIENT:  -2.6319E-01  1.6887E-01 -8.4680E-02 -4.4935E-01  5.3284E-01 -1.5811E-02  1.9143E-01  2.3768E-02 -4.0703E-02 -8.9741E-02
            -1.4521E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1727.75906570460        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2094
 NPARAMETR:  1.0324E+00  1.0005E-02  1.8421E+00  1.7274E+00  9.4045E-01  1.0430E+00  1.1725E+01  1.2966E+00  7.3922E-01  8.8795E-01
             9.4925E-01
 PARAMETER:  1.3189E-01 -4.5047E+00  7.1090E-01  6.4660E-01  3.8607E-02  1.4208E-01  2.5617E+00  3.5977E-01 -2.0216E-01 -1.8837E-02
             4.7917E-02
 GRADIENT:  -2.3730E-01  1.9178E-01  7.3925E-03 -3.9761E-01 -1.0971E-02 -1.9378E-02  2.8592E-01 -8.6568E-03 -3.1384E-03 -1.2956E-02
             4.9615E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1727.75920237812        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2275            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0325E+00  1.0000E-02  1.8421E+00  1.7274E+00  9.4046E-01  1.0430E+00  1.1717E+01  1.2967E+00  7.3922E-01  8.8667E-01
             9.4922E-01
 PARAMETER:  1.3195E-01 -4.5058E+00  7.1090E-01  6.4661E-01  3.8617E-02  1.4211E-01  2.5610E+00  3.5983E-01 -2.0216E-01 -2.0280E-02
             4.7889E-02
 GRADIENT:   6.2178E+01  2.6942E-01  1.0135E+00  1.4124E+02  8.4595E-01  6.3786E+00  3.4716E-01  3.9605E-02  2.8501E+00 -1.1562E-01
             3.0283E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1727.75936545837        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     2439
 NPARAMETR:  1.0324E+00  1.0000E-02  1.8421E+00  1.7273E+00  9.4046E-01  1.0430E+00  1.1711E+01  1.2967E+00  7.3920E-01  8.8701E-01
             9.4923E-01
 PARAMETER:  1.3191E-01 -4.5063E+00  7.1089E-01  6.4655E-01  3.8615E-02  1.4208E-01  2.5605E+00  3.5982E-01 -2.0218E-01 -1.9902E-02
             4.7893E-02
 GRADIENT:  -2.5746E-01  1.0255E-01  4.9980E-03 -4.5909E-01  1.4400E-01 -3.1526E-02  2.0346E-01 -3.1945E-02  1.9302E-02 -1.0819E-01
            -3.4947E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1727.76028450875        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2616
 NPARAMETR:  1.0325E+00  1.0000E-02  1.8417E+00  1.7273E+00  9.4039E-01  1.0430E+00  1.1572E+01  1.2964E+00  7.3905E-01  8.8721E-01
             9.4926E-01
 PARAMETER:  1.3200E-01 -4.5146E+00  7.1072E-01  6.4658E-01  3.8544E-02  1.4213E-01  2.5486E+00  3.5958E-01 -2.0240E-01 -1.9676E-02
             4.7925E-02
 GRADIENT:   8.1714E-02  0.0000E+00  1.6117E-02 -2.3793E-01  8.2694E-02  2.8884E-03  2.1018E-02 -2.0506E-03  2.0183E-02 -1.8223E-02
            -1.9350E-02

0ITERATION NO.:   83    OBJECTIVE VALUE:  -1727.76034844975        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:     2711
 NPARAMETR:  1.0325E+00  1.0000E-02  1.8416E+00  1.7272E+00  9.4030E-01  1.0430E+00  1.1462E+01  1.2964E+00  7.3900E-01  8.8731E-01
             9.4928E-01
 PARAMETER:  1.3202E-01 -4.5146E+00  7.1064E-01  6.4653E-01  3.8445E-02  1.4212E-01  2.5390E+00  3.5961E-01 -2.0246E-01 -1.9556E-02
             4.7953E-02
 GRADIENT:   7.6053E-02  0.0000E+00  5.3710E-02 -3.7001E-01 -3.0290E-02  2.6160E-02 -4.5360E-03  1.6551E-02 -1.0939E-02  1.6242E-02
            -2.5727E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2711
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1824E-04  2.2821E-04 -2.7235E-02 -5.7611E-03 -3.7286E-02
 SE:             2.9871E-02  1.9172E-03  1.7560E-02  2.9315E-02  2.0380E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9684E-01  9.0525E-01  1.2091E-01  8.4420E-01  6.7318E-02

 ETASHRINKSD(%)  1.0000E-10  9.3577E+01  4.1171E+01  1.7915E+00  3.1725E+01
 ETASHRINKVR(%)  1.0000E-10  9.9587E+01  6.5392E+01  3.5510E+00  5.3386E+01
 EBVSHRINKSD(%)  3.3841E-01  9.3984E+01  4.3732E+01  2.2152E+00  2.9359E+01
 EBVSHRINKVR(%)  6.7567E-01  9.9638E+01  6.8340E+01  4.3813E+00  5.0098E+01
 RELATIVEINF(%)  9.5994E+01  1.0501E-02  6.5202E+00  3.3932E+00  5.8093E+00
 EPSSHRINKSD(%)  4.4428E+01
 EPSSHRINKVR(%)  6.9117E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1727.7603484497545     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -992.60952188601630     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.38
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1727.760       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  1.84E+00  1.73E+00  9.40E-01  1.04E+00  1.15E+01  1.30E+00  7.39E-01  8.87E-01  9.49E-01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.13E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.68E+01  0.00E+00  4.56E+01
 
 TH 4
+       -2.37E+02  0.00E+00 -3.30E+01  7.43E+02
 
 TH 5
+       -1.53E+02  0.00E+00 -2.00E+02  3.30E+01  9.42E+02
 
 TH 6
+        1.44E+02  0.00E+00  2.20E+00 -1.03E+02 -2.95E+00  2.45E+02
 
 TH 7
+        5.90E-02  0.00E+00  6.41E-03 -6.09E-02 -2.86E-02  6.27E-02  2.10E-05
 
 TH 8
+        1.99E+01  0.00E+00  3.27E+00 -1.34E+01 -2.19E+01  1.27E+01  5.42E-03  2.58E+00
 
 TH 9
+       -1.55E+02  0.00E+00 -4.98E+00  1.01E+02  5.46E+01 -5.18E+01  1.37E-03 -8.95E+00  4.26E+02
 
 TH10
+        2.73E+00  0.00E+00  2.70E+01 -5.00E+00 -1.27E+02  7.08E+00  6.56E-03  3.35E+00  1.11E+01  1.86E+01
 
 TH11
+        1.68E+02  0.00E+00  1.97E+01 -1.28E+02 -2.03E+02  6.32E+01  4.80E-02  3.04E+01 -5.25E+01  3.30E+01  4.23E+02
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.54E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.43E+00  0.00E+00  6.49E+01
 
 TH 4
+       -6.93E+00  0.00E+00 -1.69E+01  6.53E+02
 
 TH 5
+       -2.12E+01  0.00E+00 -1.85E+02 -5.75E+01  8.41E+02
 
 TH 6
+        1.14E+01  0.00E+00  6.38E-01 -1.31E+00  1.67E+01  1.97E+02
 
 TH 7
+        1.46E-02  0.00E+00  2.09E-03 -2.87E-02 -4.19E-03  4.02E-02  6.67E-03
 
 TH 8
+        3.85E+00  0.00E+00 -1.86E+01 -2.42E+00 -1.08E+01  4.77E+00 -2.34E-02  2.51E+01
 
 TH 9
+        7.96E-01  0.00E+00  4.03E+00 -4.00E-01  1.01E+01 -7.34E+00  1.73E-02 -2.70E+00  3.40E+02
 
 TH10
+       -1.13E+01  0.00E+00  3.64E+00  2.65E+00 -1.08E+02 -3.30E+00  4.01E-02  1.65E+01  1.26E+01  8.86E+01
 
 TH11
+        6.42E+00  0.00E+00 -7.40E+00 -8.84E+00 -2.31E+01  7.29E+00  1.27E-02  1.10E+01  1.36E+01 -4.96E+00  2.46E+02
 
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
+        9.52E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.40E+00  0.00E+00  6.35E+01
 
 TH 4
+        1.17E+02  0.00E+00 -2.04E+01  6.52E+02
 
 TH 5
+        2.95E+01  0.00E+00 -1.86E+02 -8.77E+01  8.25E+02
 
 TH 6
+       -5.72E+01  0.00E+00 -4.92E+00  5.19E+01  3.60E+00  1.73E+02
 
 TH 7
+       -8.21E-04  0.00E+00  4.16E-03 -1.50E-02 -1.88E-02  1.39E-04  1.62E-05
 
 TH 8
+       -1.27E+01  0.00E+00 -1.46E+01  4.28E+00 -1.59E+01 -1.49E+00  3.34E-03  2.04E+01
 
 TH 9
+        8.49E+01  0.00E+00 -2.75E+00 -4.79E+01  1.31E+00  5.91E+00  5.08E-02  6.62E-01  2.93E+02
 
 TH10
+       -4.66E+01  0.00E+00  5.29E+00  2.47E+01 -9.45E+01  1.15E+01  1.86E-02  1.84E+01  3.72E+00  8.76E+01
 
 TH11
+       -4.01E+01  0.00E+00 -2.21E+01  2.73E+01  7.53E+01 -1.49E+01  3.62E-03  1.31E+00  3.16E+01  3.90E+00  1.60E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       41.498
Stop Time:
Sat Sep 25 13:41:04 CDT 2021
