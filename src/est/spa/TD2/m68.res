Sat Sep 25 13:42:49 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1708.55296445049        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0022E+00  1.5318E+01 -2.2314E+01  6.5025E+01  4.0151E+01  6.2449E+00  1.9385E+01  3.3073E+00  4.0637E+01  1.2667E+00
             1.9417E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1716.29400967197        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0005E+00  9.8407E-01  9.8024E-01  9.8489E-01  9.6094E-01  9.8824E-01  8.7751E-01  9.9188E-01  8.2120E-01  9.8526E-01
             9.8159E-01
 PARAMETER:  1.0047E-01  8.3945E-02  8.0039E-02  8.4774E-02  6.0157E-02  8.8166E-02 -3.0671E-02  9.1851E-02 -9.6985E-02  8.5145E-02
             8.1417E-02
 GRADIENT:   5.8774E+00  1.2338E+01 -4.8588E+00  2.1345E+01  7.4894E+00  2.1864E+00  4.6342E+00  2.2622E+00  3.1507E+00  3.4723E+00
             9.5383E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1716.91547237719        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.0146E+00  9.2058E-01  9.7440E-01  1.0277E+00  9.2216E-01  1.0132E+00  8.1404E-01  8.9083E-01  8.3735E-01  9.4126E-01
             9.6307E-01
 PARAMETER:  1.1446E-01  1.7249E-02  7.4070E-02  1.2729E-01  1.8964E-02  1.1313E-01 -1.0575E-01 -1.5601E-02 -7.7516E-02  3.9463E-02
             6.2373E-02
 GRADIENT:  -8.1507E+00  1.5966E+01  3.6676E+00  2.6033E+01 -2.2604E+00  7.4718E+00  5.4093E-01 -9.7375E-01  6.5309E+00 -3.0142E+00
             5.5860E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.47304155227        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0195E+00  8.8076E-01  8.5173E-01  1.0333E+00  8.4760E-01  9.9317E-01  9.2478E-01  7.0247E-01  7.7259E-01  8.9246E-01
             9.6279E-01
 PARAMETER:  1.1936E-01 -2.6970E-02 -6.0491E-02  1.3273E-01 -6.5345E-02  9.3142E-02  2.1796E-02 -2.5316E-01 -1.5801E-01 -1.3774E-02
             6.2077E-02
 GRADIENT:   4.2228E-01  1.8441E+00  3.2549E-01  1.4891E+00 -9.1102E-01 -4.9413E-01  5.1453E-02  3.0099E-02 -4.7940E-01  1.3515E-01
            -1.4713E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1717.55233082648        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  1.0163E+00  7.3233E-01  9.9564E-01  1.1356E+00  8.5448E-01  9.9458E-01  1.0114E+00  8.2845E-01  7.3441E-01  9.1826E-01
             9.6417E-01
 PARAMETER:  1.1619E-01 -2.1152E-01  9.5635E-02  2.2718E-01 -5.7264E-02  9.4562E-02  1.1132E-01 -8.8193E-02 -2.0869E-01  1.4729E-02
             6.3511E-02
 GRADIENT:  -1.2295E+00  7.8587E+00  6.1369E+00  9.9543E+00 -9.8140E+00  1.0473E+00  2.0629E-01 -1.1096E+00  9.8471E-01 -1.4955E-01
             1.2112E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1717.92457874900        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  1.0135E+00  5.3956E-01  1.2024E+00  1.2632E+00  8.7905E-01  9.8699E-01  1.1431E+00  1.0456E+00  6.8276E-01  9.5379E-01
             9.6526E-01
 PARAMETER:  1.1336E-01 -5.1700E-01  2.8428E-01  3.3362E-01 -2.8911E-02  8.6905E-02  2.3374E-01  1.4463E-01 -2.8161E-01  5.2692E-02
             6.4641E-02
 GRADIENT:   6.1541E-01  4.8465E+00  1.5178E+00  1.1912E+01 -3.4669E+00 -3.9961E-01 -2.3959E-01 -1.9522E-01 -3.2185E-01  2.1716E-01
             7.0455E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1718.46758413978        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      982
 NPARAMETR:  1.0093E+00  3.1143E-01  1.4339E+00  1.4162E+00  8.9165E-01  9.7909E-01  1.5999E+00  1.2798E+00  6.2044E-01  9.7446E-01
             9.6114E-01
 PARAMETER:  1.0927E-01 -1.0666E+00  4.6042E-01  4.4797E-01 -1.4685E-02  7.8872E-02  5.6997E-01  3.4668E-01 -3.7733E-01  7.4129E-02
             6.0362E-02
 GRADIENT:   8.4162E-01  5.5970E+00  5.6517E-01  2.9753E+01 -3.9134E+00 -1.7402E+00 -2.6787E-01 -1.1633E-01 -2.7638E+00  2.4987E-01
            -1.2463E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.27258725780        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1157
 NPARAMETR:  1.0046E+00  1.3832E-01  1.6928E+00  1.5264E+00  9.1495E-01  9.8286E-01  2.6368E+00  1.5596E+00  5.9277E-01  9.9207E-01
             9.6426E-01
 PARAMETER:  1.0460E-01 -1.8782E+00  6.2641E-01  5.2288E-01  1.1116E-02  8.2714E-02  1.0696E+00  5.4444E-01 -4.2294E-01  9.2037E-02
             6.3610E-02
 GRADIENT:  -2.7381E+00  1.7625E+00  1.1654E+00  1.1145E+01 -7.3553E+00  1.2699E+00  3.4028E-01  2.0407E+00  1.3824E+00  9.1318E-01
             1.2786E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1719.71010025037        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1332
 NPARAMETR:  1.0041E+00  5.9845E-02  1.7414E+00  1.5739E+00  9.1235E-01  9.7981E-01  4.0031E+00  1.5863E+00  5.7756E-01  9.8987E-01
             9.6182E-01
 PARAMETER:  1.0408E-01 -2.7160E+00  6.5470E-01  5.5355E-01  8.2717E-03  7.9606E-02  1.4871E+00  5.6140E-01 -4.4894E-01  8.9820E-02
             6.1067E-02
 GRADIENT:  -7.3666E-01  4.4957E-01 -6.8898E-01  5.4564E+00 -5.8290E-01  6.9071E-01  4.5159E-02  4.3315E-01 -1.1517E+00 -1.8907E-01
             3.5037E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1719.94244241535        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1507
 NPARAMETR:  1.0038E+00  1.0000E-02  1.8896E+00  1.6110E+00  9.3434E-01  9.7677E-01  8.1634E+00  1.7043E+00  5.7150E-01  1.0070E+00
             9.6077E-01
 PARAMETER:  1.0377E-01 -4.5190E+00  7.3636E-01  5.7683E-01  3.2082E-02  7.6495E-02  2.1997E+00  6.3318E-01 -4.5948E-01  1.0695E-01
             5.9975E-02
 GRADIENT:   2.1263E-01  0.0000E+00  3.9885E-01  7.6711E+00 -1.0578E+00 -7.9638E-02 -1.7338E-02 -3.2611E-02 -3.6601E-01 -2.2559E-01
            -2.9940E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1719.95124532645        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1683
 NPARAMETR:  1.0036E+00  1.0000E-02  1.8488E+00  1.6073E+00  9.2630E-01  9.7696E-01  8.5841E+00  1.6712E+00  5.7245E-01  1.0037E+00
             9.6120E-01
 PARAMETER:  1.0357E-01 -4.6482E+00  7.1454E-01  5.7455E-01  2.3440E-02  7.6690E-02  2.2499E+00  6.1351E-01 -4.5783E-01  1.0372E-01
             6.0430E-02
 GRADIENT:  -1.1476E-01  0.0000E+00 -2.0863E-01  1.4892E-01  1.9053E-01 -1.7418E-02 -1.5178E-02  3.6670E-02 -3.4400E-02  2.8598E-02
             1.0095E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1719.95244753298        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1859
 NPARAMETR:  1.0036E+00  1.0000E-02  1.8561E+00  1.6075E+00  9.2778E-01  9.7696E-01  9.5498E+00  1.6767E+00  5.7212E-01  1.0044E+00
             9.6118E-01
 PARAMETER:  1.0359E-01 -4.9163E+00  7.1849E-01  5.7471E-01  2.5038E-02  7.6695E-02  2.3565E+00  6.1684E-01 -4.5840E-01  1.0444E-01
             6.0409E-02
 GRADIENT:  -7.7361E-02  0.0000E+00 -3.6373E-02 -8.3095E-02  4.9697E-02 -4.4175E-02 -1.0665E-04  6.3885E-03 -2.3932E-02  2.8628E-03
            -3.4472E-03

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1719.95244947947        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1951
 NPARAMETR:  1.0036E+00  1.0000E-02  1.8572E+00  1.6076E+00  9.2799E-01  9.7702E-01  9.5512E+00  1.6775E+00  5.7216E-01  1.0046E+00
             9.6119E-01
 PARAMETER:  1.0364E-01 -4.9171E+00  7.1909E-01  5.7475E-01  2.5265E-02  7.6755E-02  2.3567E+00  6.1733E-01 -4.5834E-01  1.0455E-01
             6.0415E-02
 GRADIENT:   3.5609E-02  0.0000E+00 -4.1242E-03 -2.6496E-02  6.6876E-03 -6.1421E-03  4.4453E-05  7.7462E-04 -4.7303E-03  5.0751E-04
            -8.0002E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1951
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1001E-04  4.6863E-05 -3.5897E-02 -9.3405E-03 -4.4152E-02
 SE:             2.9869E-02  1.9811E-03  1.9034E-02  2.8837E-02  2.0217E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9706E-01  9.8113E-01  5.9300E-02  7.4601E-01  2.8972E-02

 ETASHRINKSD(%)  1.0000E-10  9.3363E+01  3.6235E+01  3.3926E+00  3.2270E+01
 ETASHRINKVR(%)  1.0000E-10  9.9560E+01  5.9341E+01  6.6702E+00  5.4126E+01
 EBVSHRINKSD(%)  3.8873E-01  9.3548E+01  3.9612E+01  3.8501E+00  2.8477E+01
 EBVSHRINKVR(%)  7.7594E-01  9.9584E+01  6.3533E+01  7.5519E+00  4.8845E+01
 RELATIVEINF(%)  9.5568E+01  9.3631E-03  9.0453E+00  2.3108E+00  8.5777E+00
 EPSSHRINKSD(%)  4.5190E+01
 EPSSHRINKVR(%)  6.9959E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.9524494794675     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.80162291572935     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.17
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.952       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.86E+00  1.61E+00  9.28E-01  9.77E-01  9.55E+00  1.68E+00  5.72E-01  1.00E+00  9.61E-01
 


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
+        1.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.35E-01  0.00E+00  4.50E+01
 
 TH 4
+       -1.53E+01  0.00E+00 -2.75E+01  1.22E+03
 
 TH 5
+       -7.07E+00  0.00E+00 -1.28E+02 -1.05E+02  6.77E+02
 
 TH 6
+       -9.03E+00  0.00E+00  5.51E-01 -2.06E+00  1.54E+01  2.14E+02
 
 TH 7
+       -6.81E-04  0.00E+00  3.84E-03  6.28E-03  1.00E-01  8.40E-02  1.88E-03
 
 TH 8
+        1.62E+00  0.00E+00 -1.57E+01 -4.83E+00 -8.53E+00 -3.50E-03 -4.35E-03  2.07E+01
 
 TH 9
+       -7.39E+00  0.00E+00  5.42E+00 -1.48E+00  1.51E+00 -3.04E+00  1.25E-01  4.10E-01  5.26E+02
 
 TH10
+       -1.08E+01  0.00E+00  1.89E+00 -3.96E+00 -1.12E+02  1.23E+01  3.63E-02  6.76E+00 -7.35E+00  6.37E+01
 
 TH11
+        7.22E+00  0.00E+00  7.67E-01 -1.79E+01 -1.42E+01  9.16E+00  3.36E-02  4.32E+00  2.03E+01  4.37E+01  2.41E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.097
Stop Time:
Sat Sep 25 13:43:21 CDT 2021
