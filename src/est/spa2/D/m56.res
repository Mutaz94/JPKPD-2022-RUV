Thu Sep 30 09:19:10 CDT 2021
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
$DATA ../../../../data/spa2/D/dat56.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m56.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26081.6214906634        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0661E+02  4.8268E+02 -1.6806E+01  2.6047E+02  1.9515E+02 -2.5434E+03 -9.4163E+02 -7.1127E+01 -1.6843E+03 -7.4863E+02
            -5.0589E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -568.317194200726        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.5879E+00  1.2163E+00  9.3311E-01  1.5500E+00  1.0556E+00  2.2729E+00  1.5165E+00  9.8804E-01  1.4140E+00  1.1020E+00
             1.3923E+01
 PARAMETER:  5.6242E-01  2.9583E-01  3.0769E-02  5.3827E-01  1.5414E-01  9.2104E-01  5.1639E-01  8.7969E-02  4.4639E-01  1.9710E-01
             2.7336E+00
 GRADIENT:   6.8962E+01 -1.2708E+01 -8.1720E+00  6.2983E+00  2.1557E+01  4.2082E+01 -2.2537E+01  2.5665E+00 -3.1884E+01  1.2908E+01
             1.2461E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -645.708706721282        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4779E+00  2.8620E+00  3.8526E+00  1.0156E+00  3.2850E+00  3.0571E+00  5.8722E+00  4.8547E-01  4.2636E+00  5.7490E-01
             1.2546E+01
 PARAMETER:  4.9059E-01  1.1515E+00  1.4487E+00  1.1544E-01  1.2894E+00  1.2175E+00  1.8702E+00 -6.2264E-01  1.5501E+00 -4.5357E-01
             2.6294E+00
 GRADIENT:   3.2387E+01  2.3348E+01 -7.9339E+00 -1.0794E+01  1.0710E+01  8.1662E+01  1.1939E+02  8.6048E-03  3.7124E+01  8.5904E-01
             1.9866E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -699.827419447030        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.3106E+00  8.3841E-01  7.8458E+00  1.7979E+00  1.6829E+00  2.4634E+00  5.5383E+00  6.0157E-01  1.7045E+00  9.7994E-01
             1.2340E+01
 PARAMETER:  3.7052E-01 -7.6248E-02  2.1600E+00  6.8661E-01  6.2051E-01  1.0015E+00  1.8117E+00 -4.0821E-01  6.3324E-01  7.9739E-02
             2.6128E+00
 GRADIENT:   1.4771E+01 -3.9671E+00  3.5906E+00  1.9635E+01 -4.2011E+01  5.7765E+01  1.7089E+01  8.1781E-02  1.3439E+01  7.6219E+00
             1.8316E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -731.184378765193        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.2310E+00  9.1818E-01  1.5091E+01  1.8715E+00  2.6584E+00  1.9762E+00  6.7232E+00  1.5394E+00  1.8082E+00  1.0028E+00
             1.2040E+01
 PARAMETER:  3.0781E-01  1.4643E-02  2.8141E+00  7.2673E-01  1.0777E+00  7.8116E-01  2.0056E+00  5.3142E-01  6.9235E-01  1.0281E-01
             2.5882E+00
 GRADIENT:  -2.2958E+00  3.2848E+00 -1.0137E+00 -5.6326E+00 -2.5436E+00 -2.8374E+00 -2.0513E+00  8.1431E-02  2.1754E+01  4.4834E+00
             1.2766E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -746.745792277379        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  1.1804E+00  8.0140E-01  1.7589E+01  1.4660E+00  3.1215E+00  2.0095E+00  6.8435E+00  6.1626E+00  8.7583E-01  8.2964E-01
             1.0556E+01
 PARAMETER:  2.6582E-01 -1.2140E-01  2.9673E+00  4.8257E-01  1.2383E+00  7.9788E-01  2.0233E+00  1.9185E+00 -3.2581E-02 -8.6758E-02
             2.4567E+00
 GRADIENT:   6.2312E+00 -7.1908E+00 -1.0374E+00 -1.8710E+01  6.8821E+00  1.7706E+00 -9.4298E-01  4.0474E-01  4.6140E+00  2.3679E+00
             1.2994E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -749.026948638587        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  1.1539E+00  1.0472E+00  1.3780E+01  1.4426E+00  2.9041E+00  1.9870E+00  6.8417E+00  1.6262E+00  7.2850E-01  3.5843E-01
             1.0602E+01
 PARAMETER:  2.4312E-01  1.4611E-01  2.7232E+00  4.6645E-01  1.1661E+00  7.8661E-01  2.0230E+00  5.8626E-01 -2.1677E-01 -9.2603E-01
             2.4611E+00
 GRADIENT:   5.4833E+00  3.3929E+00 -7.7508E-01  1.7701E+01  5.1370E+00  1.4847E+01  9.1221E+01  6.0667E-02 -7.8587E-01  4.6815E-01
             3.6753E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -749.388344786768        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      850
 NPARAMETR:  1.1648E+00  9.8647E-01  1.3805E+01  1.4423E+00  2.9043E+00  1.9956E+00  6.8059E+00  6.0940E-01  7.6631E-01  1.3614E-01
             1.0581E+01
 PARAMETER:  2.5252E-01  8.6373E-02  2.7250E+00  4.6621E-01  1.1662E+00  7.9093E-01  2.0178E+00 -3.9529E-01 -1.6617E-01 -1.8941E+00
             2.4591E+00
 GRADIENT:  -9.2786E-01 -4.5260E-02 -7.9271E-01  1.0762E+00  5.8924E+00 -6.5904E-01  4.4043E+00  8.4760E-03  4.9402E-01  6.6724E-02
             8.1302E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -749.862128544312        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1029             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1657E+00  9.9071E-01  2.3820E+01  1.4330E+00  2.6173E+00  1.9998E+00  6.8516E+00  1.9442E-01  7.6120E-01  3.1692E-02
             1.0536E+01
 PARAMETER:  2.5331E-01  9.0670E-02  3.2705E+00  4.5979E-01  1.0622E+00  7.9305E-01  2.0245E+00 -1.5377E+00 -1.7286E-01 -3.3517E+00
             2.4548E+00
 GRADIENT:   1.2814E+01  1.6162E+00 -4.4862E-02  1.3147E+01 -2.5048E+00  1.9242E+01  9.4071E+01  3.0490E-04 -4.6858E-01  4.7723E-03
             3.2344E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -749.982338711909        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1208
 NPARAMETR:  1.1643E+00  9.7927E-01  4.5800E+01  1.4343E+00  2.7140E+00  2.0004E+00  6.8247E+00  1.6224E-01  7.7503E-01  1.2027E-02
             1.0499E+01
 PARAMETER:  2.5209E-01  7.9052E-02  3.9243E+00  4.6070E-01  1.0984E+00  7.9333E-01  2.0205E+00 -1.7187E+00 -1.5485E-01 -4.3206E+00
             2.4513E+00
 GRADIENT:   1.6501E-01  2.3245E-01 -3.2055E-02  1.1827E+00 -1.1178E+00  1.0897E+00  4.8667E+00  5.0277E-05 -8.1067E-02  6.3647E-04
             8.4413E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -750.030559112891        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1391
 NPARAMETR:  1.1649E+00  9.7308E-01  1.9026E+02  1.4377E+00  2.7673E+00  2.0007E+00  6.9111E+00  1.4419E-01  7.7138E-01  1.0000E-02
             1.0491E+01
 PARAMETER:  2.5262E-01  7.2714E-02  5.3484E+00  4.6307E-01  1.1179E+00  7.9352E-01  2.0331E+00 -1.8366E+00 -1.5958E-01 -4.6941E+00
             2.4505E+00
 GRADIENT:   4.1348E-01  4.0855E-01 -4.1425E-03 -1.9056E-01 -5.3123E-01  9.0443E-01  7.2135E+00  1.9571E-06 -1.2363E-01  0.0000E+00
             5.6461E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -750.041553512056        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1570
 NPARAMETR:  1.1646E+00  9.6311E-01  1.0644E+03  1.4424E+00  2.7864E+00  1.9994E+00  6.9313E+00  1.4474E-01  7.8281E-01  1.0000E-02
             1.0486E+01
 PARAMETER:  2.5237E-01  6.2414E-02  7.0701E+00  4.6629E-01  1.1247E+00  7.9285E-01  2.0360E+00 -1.8328E+00 -1.4487E-01 -4.6941E+00
             2.4500E+00
 GRADIENT:   3.5308E-01  2.4618E-01 -8.0582E-04 -1.0713E+00 -1.0691E-01  5.9410E-01  7.4434E+00  5.1385E-08  1.3268E-01  0.0000E+00
             4.7609E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -750.046612409818        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1753             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1634E+00  9.5000E-01  4.9075E+03  1.4480E+00  2.7931E+00  1.9993E+00  6.9501E+00  1.4381E-01  7.8390E-01  1.0000E-02
             1.0475E+01
 PARAMETER:  2.5136E-01  4.8709E-02  8.5985E+00  4.7020E-01  1.1271E+00  7.9279E-01  2.0387E+00 -1.8393E+00 -1.4347E-01 -4.6941E+00
             2.4490E+00
 GRADIENT:   1.2443E+01  7.6663E-01 -1.8467E-04  9.8790E+00  6.3637E-01  1.8654E+01  9.9103E+01  3.2581E-07 -2.9606E-02  0.0000E+00
             2.8390E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -750.047608344270        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1936             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1641E+00  9.4748E-01  5.8997E+03  1.4498E+00  2.7935E+00  1.9999E+00  6.9573E+00  2.1645E-01  7.8762E-01  1.0000E-02
             1.0481E+01
 PARAMETER:  2.5191E-01  4.6055E-02  8.7826E+00  4.7141E-01  1.1273E+00  7.9308E-01  2.0398E+00 -1.4304E+00 -1.3874E-01 -4.6941E+00
             2.4495E+00
 GRADIENT:   1.2648E+01  7.2923E-01 -1.5886E-04  9.6666E+00  6.7065E-01  1.8778E+01  9.9328E+01  1.1027E-07  5.0751E-02  0.0000E+00
             2.9038E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -750.047949696086        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2123
 NPARAMETR:  1.1638E+00  9.4627E-01  1.8490E+04  1.4505E+00  2.7933E+00  1.9999E+00  6.9614E+00  2.1810E-01  7.8833E-01  1.0000E-02
             1.0480E+01
 PARAMETER:  2.5171E-01  4.4771E-02  9.9250E+00  4.7193E-01  1.1272E+00  7.9309E-01  2.0404E+00 -1.4228E+00 -1.3783E-01 -4.6941E+00
             2.4495E+00
 GRADIENT:  -4.2591E-02  5.2335E-02 -5.2296E-05  1.5348E-02 -4.9018E-02  9.8114E-01  7.2148E+00  3.7372E-07 -9.8918E-02  0.0000E+00
            -3.9552E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -750.048327333017        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2317             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1638E+00  9.4147E-01  1.8754E+04  1.4514E+00  2.7946E+00  1.9999E+00  6.9680E+00  2.2564E-01  7.9136E-01  1.0000E-02
             1.0481E+01
 PARAMETER:  2.5167E-01  3.9686E-02  9.9392E+00  4.7255E-01  1.1277E+00  7.9309E-01  2.0413E+00 -1.3888E+00 -1.3401E-01 -4.6941E+00
             2.4496E+00
 GRADIENT:   1.2267E+01  5.8806E-01 -5.1878E-05  9.3338E+00  7.4992E-01  1.8992E+01  9.9807E+01  7.4101E-07  1.2503E-01  0.0000E+00
             2.9433E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -750.048528942987        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2503
 NPARAMETR:  1.1638E+00  9.4142E-01  1.9233E+04  1.4523E+00  2.7936E+00  1.9998E+00  6.9718E+00  2.2066E-01  7.9373E-01  1.0000E-02
             1.0480E+01
 PARAMETER:  2.5168E-01  3.9639E-02  9.9644E+00  4.7317E-01  1.1273E+00  7.9304E-01  2.0419E+00 -1.4111E+00 -1.3101E-01 -4.6941E+00
             2.4495E+00
 GRADIENT:  -1.2221E-01  1.4557E-02 -5.3653E-05 -4.2018E-01  2.5966E-02  7.8463E-01  7.3833E+00  7.9750E-08  4.0470E-02  0.0000E+00
            -1.4037E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -750.048710255364        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2688
 NPARAMETR:  1.1639E+00  9.4067E-01  2.2974E+04  1.4554E+00  2.7911E+00  1.9999E+00  6.9757E+00  2.2541E-01  7.9347E-01  1.0000E-02
             1.0478E+01
 PARAMETER:  2.5170E-01  3.8742E-02  1.0225E+01  4.7364E-01  1.1274E+00  7.9304E-01  2.0421E+00 -1.3898E+00 -1.3004E-01 -4.6941E+00
             2.4495E+00
 GRADIENT:  -1.3700E-02 -3.5944E-03  1.2749E-04 -3.3088E-01  2.6645E-02 -4.5954E-03 -2.5888E-02  8.5921E-07  4.7092E-02  0.0000E+00
             1.1052E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2688
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0110E-02  3.6157E-02  6.9887E-08 -6.4058E-02  6.3621E-06
 SE:             2.8213E-02  2.3690E-02  2.6938E-08  1.1497E-02  3.9881E-05
 N:                     100         100         100         100         100

 P VAL.:         4.7598E-01  1.2695E-01  9.4767E-03  2.5297E-08  8.7326E-01

 ETASHRINKSD(%)  5.4832E+00  2.0635E+01  1.0000E+02  6.1484E+01  9.9866E+01
 ETASHRINKVR(%)  1.0666E+01  3.7012E+01  1.0000E+02  8.5165E+01  1.0000E+02
 EBVSHRINKSD(%)  7.5488E+00  1.4676E+01  1.0000E+02  6.5975E+01  9.9802E+01
 EBVSHRINKVR(%)  1.4528E+01  2.7198E+01  1.0000E+02  8.8423E+01  1.0000E+02
 RELATIVEINF(%)  8.2892E+01  4.1071E+01  0.0000E+00  5.7281E+00  6.2477E-05
 EPSSHRINKSD(%)  6.9168E+00
 EPSSHRINKVR(%)  1.3355E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -750.04871025536386     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       352.67752959024324     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    67.28
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -750.049       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  9.41E-01  2.50E+04  1.45E+00  2.79E+00  2.00E+00  6.97E+00  2.25E-01  7.95E-01  1.00E-02  1.05E+01
 


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
+        2.52E+02
 
 TH 2
+       -4.82E+00  2.96E+00
 
 TH 3
+       -1.02E-06  1.29E-07  2.26E-14
 
 TH 4
+       -1.61E+02  2.05E+01  2.10E-06  2.52E+02
 
 TH 5
+        9.45E+00 -1.10E+00 -1.37E-07 -1.51E+01  9.30E-01
 
 TH 6
+       -6.65E+01 -2.31E+00  8.88E-07  6.25E+01 -4.74E+00  6.26E+01
 
 TH 7
+        9.62E+00 -9.78E-01 -1.15E-07 -1.35E+01  8.18E-01 -3.98E+00  7.33E-01
 
 TH 8
+        2.08E+00 -5.07E-02 -1.23E-08 -1.59E+00  9.72E-02 -7.25E-01  9.38E-02  1.81E-02
 
 TH 9
+        6.03E+01 -6.90E+00 -7.77E-07 -9.08E+01  5.50E+00 -2.55E+01  4.90E+00  5.98E-01  3.29E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.50E-01 -1.29E+00 -6.83E-08 -9.04E+00  5.10E-01  2.07E-01  4.23E-01  8.20E-03  3.06E+00  0.00E+00  9.56E-01
 
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
+        1.67E+02
 
 TH 2
+       -3.90E+00  3.14E+01
 
 TH 3
+        1.58E-07  3.76E-06  1.11E-12
 
 TH 4
+       -8.66E+00  2.20E+01  2.43E-07  1.46E+02
 
 TH 5
+       -5.05E-01 -1.44E+00  6.21E-08 -7.91E+00  5.55E+00
 
 TH 6
+        1.50E-01 -1.63E+00  3.41E-07  3.49E+00 -1.07E+00  4.42E+01
 
 TH 7
+        4.89E-01  2.37E+00  2.30E-08 -7.20E+00  2.10E-01 -4.39E-01  2.21E+00
 
 TH 8
+        5.83E-01  2.37E+00  1.18E-06 -3.79E-01 -3.51E-02  2.15E-02  3.00E-02  1.05E+00
 
 TH 9
+       -7.19E-01  3.78E+00 -2.72E-06 -4.62E+01  2.78E+00 -9.59E-01  1.81E+00 -3.07E-02  2.75E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.28E+00 -1.98E+00 -9.36E-09 -9.75E+00  3.95E-01  4.10E-01  3.33E-01 -9.19E-03  3.09E+00  0.00E+00  5.95E+00
 
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
+        1.79E+02
 
 TH 2
+        3.56E+01  1.95E+01
 
 TH 3
+       -1.89E-09  1.02E-09  5.83E-17
 
 TH 4
+        7.86E+01  2.37E+01  1.77E-08  1.45E+02
 
 TH 5
+       -6.00E+00 -2.66E-01 -8.56E-09 -1.17E+01  2.58E+00
 
 TH 6
+        3.95E+01  5.58E+00 -2.54E-09 -1.30E+01  1.87E-01  4.79E+01
 
 TH 7
+       -2.99E+00  3.17E+00 -1.01E-09 -8.19E+00  1.51E+00  1.43E+00  3.13E+00
 
 TH 8
+        4.87E-03  9.00E-04  1.37E-12 -1.15E-03 -1.87E-04 -3.36E-03  1.09E-04  1.39E-05
 
 TH 9
+       -8.98E+00 -1.66E+00 -5.26E-09 -3.62E+01  2.92E+00  4.46E+00  1.44E+00  1.78E-03  1.96E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.05E+01 -3.59E+00 -4.85E-09 -1.62E+01  2.24E+00 -9.94E+00  2.07E+00  3.82E-03  4.13E+00  0.00E+00  8.53E+01
 
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
 #CPUT: Total CPU Time in Seconds,       79.211
Stop Time:
Thu Sep 30 09:20:31 CDT 2021
