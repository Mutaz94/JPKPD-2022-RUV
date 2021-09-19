Sat Sep 18 12:26:37 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1663.22333416684        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0190E+01 -5.0047E+01 -1.1514E+01 -6.4353E+01  2.0472E+01  2.2804E+01 -6.2455E+00  7.3569E+00 -1.0744E+01  8.1262E+00
            -9.0917E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1666.40169086174        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.9633E-01  1.0234E+00  1.0140E+00  1.0325E+00  9.9250E-01  9.3485E-01  1.0203E+00  9.7086E-01  1.0297E+00  9.6404E-01
             1.0291E+00
 PARAMETER:  9.6318E-02  1.2312E-01  1.1393E-01  1.3202E-01  9.2468E-02  3.2627E-02  1.2009E-01  7.0426E-02  1.2926E-01  6.3381E-02
             1.2866E-01
 GRADIENT:  -5.9992E-01 -1.1027E+00  8.8166E-01 -3.6027E+00  9.2149E-01 -6.4327E+00 -2.6118E+00  4.6232E+00 -1.1308E+00  3.5384E+00
             1.1856E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1668.71906975322        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  9.9577E-01  8.9242E-01  8.3520E-01  1.1098E+00  8.3783E-01  9.5189E-01  1.2043E+00  6.1792E-01  9.5111E-01  7.7664E-01
             1.0260E+00
 PARAMETER:  9.5763E-02 -1.3815E-02 -8.0080E-02  2.0417E-01 -7.6941E-02  5.0689E-02  2.8593E-01 -3.8139E-01  4.9879E-02 -1.5277E-01
             1.2567E-01
 GRADIENT:  -4.1724E+00  6.7295E+00 -3.2390E+00  1.6770E+01  9.1117E+00  2.4779E-01 -2.7115E+00  1.2309E+00 -2.3272E+00 -3.8282E+00
            -3.0626E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1669.39170365900        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      489
 NPARAMETR:  9.9760E-01  7.8442E-01  7.9225E-01  1.1627E+00  7.6579E-01  9.5145E-01  1.3574E+00  4.2733E-01  9.1317E-01  7.6809E-01
             1.0260E+00
 PARAMETER:  9.7593E-02 -1.4281E-01 -1.3288E-01  2.5074E-01 -1.6685E-01  5.0227E-02  4.0559E-01 -7.5021E-01  9.1654E-03 -1.6385E-01
             1.2562E-01
 GRADIENT:   9.5031E-01  7.9624E+00  6.7065E+00  5.4112E+00 -1.0498E+01  1.5849E-01 -1.7079E-01 -2.3528E-01 -5.8302E-01 -6.3867E-01
             9.4618E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1669.80086927951        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      666
 NPARAMETR:  9.9483E-01  5.6558E-01  8.4554E-01  1.2926E+00  7.2158E-01  9.4705E-01  1.7288E+00  4.7959E-01  8.4739E-01  7.8057E-01
             1.0292E+00
 PARAMETER:  9.4817E-02 -4.6991E-01 -6.7778E-02  3.5669E-01 -2.2632E-01  4.5594E-02  6.4744E-01 -6.3483E-01 -6.5597E-02 -1.4773E-01
             1.2875E-01
 GRADIENT:   1.0494E+00  5.2097E+00  5.2208E+00  4.3236E+00 -8.5443E+00 -5.5431E-01  1.2852E+00 -3.6257E-01 -1.9333E+00  6.2043E-01
             1.5605E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1670.11320335816        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      841
 NPARAMETR:  9.9012E-01  3.5942E-01  9.0065E-01  1.4157E+00  6.9452E-01  9.4489E-01  2.2945E+00  5.6985E-01  8.0689E-01  7.9028E-01
             1.0269E+00
 PARAMETER:  9.0074E-02 -9.2327E-01 -4.6388E-03  4.4766E-01 -2.6453E-01  4.3316E-02  9.3051E-01 -4.6239E-01 -1.1457E-01 -1.3536E-01
             1.2657E-01
 GRADIENT:  -3.6872E-01  2.7612E+00  4.2617E+00  9.9524E+00 -6.0092E+00 -9.9864E-03  2.8790E-01 -2.9199E-01 -1.0985E+00  1.8349E-02
             5.4807E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1670.19753025210        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  9.8802E-01  2.6509E-01  9.0623E-01  1.4638E+00  6.7595E-01  9.4338E-01  2.7129E+00  5.9701E-01  7.9257E-01  7.9228E-01
             1.0252E+00
 PARAMETER:  8.7947E-02 -1.2277E+00  1.5360E-03  4.8102E-01 -2.9163E-01  4.1718E-02  1.0980E+00 -4.1583E-01 -1.3248E-01 -1.3285E-01
             1.2488E-01
 GRADIENT:  -2.8343E-01  6.3929E-01  8.6853E-01  2.7345E+00 -2.3415E+00  9.4180E-02  8.3084E-02 -8.4564E-03 -6.6512E-01  3.3215E-01
             1.8963E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1670.25077897865        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  9.8520E-01  1.6423E-01  9.5490E-01  1.5266E+00  6.7651E-01  9.4092E-01  3.4902E+00  6.8131E-01  7.7713E-01  7.9443E-01
             1.0241E+00
 PARAMETER:  8.5094E-02 -1.7065E+00  5.3850E-02  5.2302E-01 -2.9080E-01  3.9104E-02  1.3500E+00 -2.8374E-01 -1.5215E-01 -1.3013E-01
             1.2377E-01
 GRADIENT:  -4.0384E-01  9.5473E-01  2.7639E+00  6.6281E+00 -3.5132E+00  9.8951E-03  3.6442E-01 -2.6330E-01 -6.3780E-01 -5.0726E-01
            -4.4513E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1670.71940592814        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1371
 NPARAMETR:  9.8199E-01  5.5382E-02  9.6012E-01  1.5782E+00  6.5671E-01  9.4191E-01  5.7452E+00  7.3110E-01  7.6047E-01  7.7088E-01
             1.0261E+00
 PARAMETER:  8.1825E-02 -2.7935E+00  5.9305E-02  5.5630E-01 -3.2052E-01  4.0159E-02  1.8484E+00 -2.1321E-01 -1.7382E-01 -1.6023E-01
             1.2578E-01
 GRADIENT:  -1.1201E+00  5.8131E+00 -6.0583E+00 -1.9128E+01  8.2865E+00  6.4657E-01  1.2449E+01 -9.3965E-01 -1.0029E+01 -2.7252E+00
             1.8566E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1671.43240188040        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1546
 NPARAMETR:  9.8152E-01  1.0159E-02  9.1083E-01  1.5897E+00  6.2700E-01  9.3826E-01  1.1584E+01  6.7499E-01  7.5795E-01  7.4203E-01
             1.0250E+00
 PARAMETER:  8.1350E-02 -4.4894E+00  6.6045E-03  5.6354E-01 -3.6681E-01  3.6271E-02  2.5496E+00 -2.9306E-01 -1.7714E-01 -1.9836E-01
             1.2472E-01
 GRADIENT:  -1.3703E-01 -2.1896E+00  3.1874E+00  9.2873E+00 -3.9561E+00  1.1673E-01 -5.0897E+00 -1.1209E-02  3.9101E+00  4.6038E-01
            -1.5338E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1671.78248337820        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1728
 NPARAMETR:  9.8145E-01  1.0980E-02  8.8932E-01  1.5811E+00  6.1747E-01  9.3872E-01  1.1993E+01  6.6194E-01  7.5391E-01  7.4287E-01
             1.0264E+00
 PARAMETER:  8.1278E-02 -4.4117E+00 -1.7294E-02  5.5811E-01 -3.8212E-01  3.6758E-02  2.5843E+00 -3.1259E-01 -1.8249E-01 -1.9724E-01
             1.2607E-01
 GRADIENT:  -2.6186E-01  3.4023E-01 -2.0491E+00 -5.7687E+00  8.2289E-01  1.4588E-01  7.6105E-01  5.6574E-01 -1.2516E+00  9.8918E-01
             8.1735E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1671.79624166606        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1904
 NPARAMETR:  9.8155E-01  1.1418E-02  9.0613E-01  1.5862E+00  6.2535E-01  9.3825E-01  1.1838E+01  6.7411E-01  7.5511E-01  7.4130E-01
             1.0253E+00
 PARAMETER:  8.1380E-02 -4.3725E+00  1.4297E-03  5.6134E-01 -3.6945E-01  3.6257E-02  2.5713E+00 -2.9436E-01 -1.8090E-01 -1.9935E-01
             1.2500E-01
 GRADIENT:   1.7589E-02  1.1840E-01 -2.7150E-02  8.1888E-02  2.1332E-01 -3.0930E-02  2.3574E-01 -8.0460E-02 -1.0318E-01 -1.4711E-01
            -7.5649E-02

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1671.79664188035        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1996
 NPARAMETR:  9.8155E-01  1.1386E-02  9.0503E-01  1.5859E+00  6.2481E-01  9.3830E-01  1.1851E+01  6.7358E-01  7.5507E-01  7.4137E-01
             1.0254E+00
 PARAMETER:  8.1382E-02 -4.3754E+00  2.1403E-04  5.6114E-01 -3.7030E-01  3.6311E-02  2.5724E+00 -2.9515E-01 -1.8094E-01 -1.9925E-01
             1.2506E-01
 GRADIENT:  -1.7390E-03 -2.3720E-03  5.0114E-03  4.2560E-02 -1.3030E-02 -7.5197E-03 -6.0329E-03 -1.1841E-03  1.2465E-02  2.1858E-03
            -2.5356E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1996
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2620E-04  1.4982E-02 -1.5995E-02 -7.6295E-03 -1.7030E-02
 SE:             2.9809E-02  7.2479E-03  1.5228E-02  2.8784E-02  2.1916E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9395E-01  3.8730E-02  2.9357E-01  7.9097E-01  4.3713E-01

 ETASHRINKSD(%)  1.3641E-01  7.5719E+01  4.8983E+01  3.5683E+00  2.6578E+01
 ETASHRINKVR(%)  2.7262E-01  9.4104E+01  7.3972E+01  7.0093E+00  4.6092E+01
 EBVSHRINKSD(%)  4.8527E-01  8.0996E+01  4.9964E+01  3.5119E+00  2.4663E+01
 EBVSHRINKVR(%)  9.6819E-01  9.6388E+01  7.4964E+01  6.9005E+00  4.3243E+01
 RELATIVEINF(%)  9.8951E+01  2.9790E+00  2.4017E+00  6.5832E+01  5.4477E+00
 EPSSHRINKSD(%)  4.4081E+01
 EPSSHRINKVR(%)  6.8731E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1671.7966418803480     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -936.64581531660986     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.40
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1671.797       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.82E-01  1.14E-02  9.05E-01  1.59E+00  6.25E-01  9.38E-01  1.19E+01  6.74E-01  7.55E-01  7.41E-01  1.03E+00
 


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
+        1.28E+03
 
 TH 2
+       -8.27E+02  1.03E+08
 
 TH 3
+       -1.11E+02  1.05E+04  2.52E+03
 
 TH 4
+        4.07E+01  1.69E+04 -7.84E+02  3.23E+05
 
 TH 5
+       -2.10E+02 -2.57E+04  1.71E+03  1.59E+03  4.78E+06
 
 TH 6
+       -5.85E+00  1.01E+03  1.59E+02 -6.83E+01  2.80E+02  2.45E+02
 
 TH 7
+       -1.22E+00 -6.21E+02  1.64E+01  2.66E+01 -4.05E+01  1.59E+00  2.75E+02
 
 TH 8
+       -3.09E+01  1.89E+03  3.84E+02 -1.40E+02  6.42E+02  1.65E+01  3.00E+00  1.29E+02
 
 TH 9
+       -2.07E+02  1.29E+04  3.64E+03 -9.11E+02 -8.09E+06  2.42E+02  2.06E+01  8.17E+02  1.37E+07
 
 TH10
+       -4.51E+01  2.92E+03  8.83E+02 -2.43E+02  1.01E+03  7.10E+01  4.60E+00  2.65E+02  1.85E+03  6.26E+02
 
 TH11
+       -4.80E+01  5.06E+02  1.34E+02 -4.67E+01  1.81E+02  3.54E+01  8.28E-01  6.23E+01  2.59E+02  9.94E+01  2.15E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.260
Stop Time:
Sat Sep 18 12:27:11 CDT 2021
