Sat Sep 25 10:04:34 CDT 2021
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
$DATA ../../../../data/spa/S1/dat74.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1676.45084662625        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0847E+01 -4.7988E+01 -1.3354E+01 -6.0913E+01  1.0390E+01  2.3967E+01 -6.1200E+00  8.6630E+00 -8.3503E+00  1.3804E+01
             6.0390E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1679.92641195773        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9270E-01  1.0269E+00  1.0327E+00  1.0297E+00  1.0101E+00  9.4918E-01  1.0425E+00  9.4335E-01  1.0387E+00  9.2398E-01
             9.9761E-01
 PARAMETER:  9.2671E-02  1.2659E-01  1.3222E-01  1.2929E-01  1.1005E-01  4.7841E-02  1.4166E-01  4.1679E-02  1.3795E-01  2.0935E-02
             9.7610E-02
 GRADIENT:  -8.0133E+00 -3.1760E-01  2.2901E+00 -7.2289E-01  5.7834E+00  7.4919E-01 -1.7453E+00  2.8769E+00  2.2138E+00 -5.0436E-01
             1.4654E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1681.68559840459        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      332
 NPARAMETR:  9.9457E-01  8.8424E-01  8.3618E-01  1.1163E+00  8.3981E-01  9.4541E-01  1.2416E+00  6.1303E-01  9.3544E-01  7.5921E-01
             1.0036E+00
 PARAMETER:  9.4551E-02 -2.3027E-02 -7.8917E-02  2.1004E-01 -7.4577E-02  4.3862E-02  3.1636E-01 -3.8933E-01  3.3257E-02 -1.7548E-01
             1.0361E-01
 GRADIENT:  -6.1105E+00  8.8860E+00 -8.8259E+00  2.4609E+01  1.1111E+01 -1.2030E+00 -1.2352E+00  1.4313E+00 -2.3713E+00 -1.6280E+00
             4.7753E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1682.24141110953        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.9699E-01  7.6760E-01  8.2312E-01  1.1658E+00  7.8321E-01  9.4647E-01  1.3911E+00  4.9015E-01  9.0598E-01  7.4604E-01
             9.9198E-01
 PARAMETER:  9.6987E-02 -1.6448E-01 -9.4651E-02  2.5339E-01 -1.4435E-01  4.4987E-02  4.3008E-01 -6.1305E-01  1.2610E-03 -1.9297E-01
             9.1946E-02
 GRADIENT:   2.0336E+00  1.8150E+00  4.0958E+00 -4.7321E+00 -5.6827E+00 -5.2447E-01  1.7355E-01 -4.7370E-01  2.2253E-01 -1.1183E+00
            -2.2389E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1682.67830295530        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      685
 NPARAMETR:  9.9269E-01  5.5032E-01  9.1645E-01  1.3064E+00  7.6089E-01  9.4587E-01  1.7459E+00  6.1410E-01  8.4276E-01  7.6599E-01
             9.8857E-01
 PARAMETER:  9.2664E-02 -4.9726E-01  1.2748E-02  3.6726E-01 -1.7326E-01  4.4346E-02  6.5729E-01 -3.8760E-01 -7.1074E-02 -1.6658E-01
             8.8508E-02
 GRADIENT:  -3.6889E-01  2.3957E+00 -5.9647E-01  5.6248E+00  6.0427E-01  4.3315E-01 -2.2821E-02  2.3526E-01 -1.3707E+00 -3.2150E-01
            -1.2202E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1682.89642606268        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      863
 NPARAMETR:  9.8990E-01  3.9108E-01  9.3816E-01  1.3984E+00  7.2404E-01  9.4181E-01  2.1734E+00  6.2003E-01  8.1357E-01  7.7755E-01
             9.9271E-01
 PARAMETER:  8.9848E-02 -8.3885E-01  3.6168E-02  4.3529E-01 -2.2291E-01  4.0051E-02  8.7627E-01 -3.7798E-01 -1.0633E-01 -1.5161E-01
             9.2685E-02
 GRADIENT:  -2.4956E-01  2.6207E+00  5.1550E+00  5.4388E+00 -8.8525E+00 -2.4537E-01  1.6935E-01 -1.3001E-01  3.4851E-01  7.5289E-01
             5.0240E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1683.05545253524        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1039
 NPARAMETR:  9.8671E-01  2.5190E-01  9.8692E-01  1.4779E+00  7.1595E-01  9.3911E-01  2.7785E+00  6.8140E-01  7.8971E-01  7.8639E-01
             9.9323E-01
 PARAMETER:  8.6624E-02 -1.2787E+00  8.6835E-02  4.9062E-01 -2.3414E-01  3.7172E-02  1.1219E+00 -2.8361E-01 -1.3609E-01 -1.4030E-01
             9.3208E-02
 GRADIENT:  -7.2596E-02  3.5238E-01  2.3336E+00 -2.9365E-01 -1.1796E+00 -2.6717E-01  1.3104E-01 -7.0881E-01  6.4620E-01 -6.2178E-01
             3.7881E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1683.18026744736        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1214
 NPARAMETR:  9.8420E-01  1.5015E-01  1.0247E+00  1.5403E+00  7.0964E-01  9.3785E-01  3.6128E+00  7.6111E-01  7.7115E-01  7.9067E-01
             9.9031E-01
 PARAMETER:  8.4069E-02 -1.7961E+00  1.2440E-01  5.3200E-01 -2.4299E-01  3.5831E-02  1.3845E+00 -1.7298E-01 -1.5987E-01 -1.3487E-01
             9.0259E-02
 GRADIENT:  -7.2118E-02  4.6532E-01  3.7689E-01  4.2408E+00 -1.3325E+00  8.0815E-02 -6.0423E-02  7.9574E-02 -3.1061E-01  1.0396E-01
            -2.8364E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1683.75151065131        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1393
 NPARAMETR:  9.8054E-01  2.2521E-02  1.0196E+00  1.6063E+00  6.7600E-01  9.3420E-01  8.5900E+00  8.5848E-01  7.5890E-01  7.3584E-01
             9.9270E-01
 PARAMETER:  8.0347E-02 -3.6933E+00  1.1941E-01  5.7393E-01 -2.9156E-01  3.1937E-02  2.2506E+00 -5.2587E-02 -1.7589E-01 -2.0674E-01
             9.2676E-02
 GRADIENT:  -1.8317E+00 -2.3796E+00  5.9075E+00  2.3532E+01 -1.7911E+01 -4.0612E-01 -6.5662E+00  4.3876E+00  6.4212E+00  3.5811E+00
             2.1411E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1684.40525236859        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  9.8098E-01  1.2172E-02  9.6135E-01  1.5935E+00  6.5391E-01  9.3543E-01  1.1439E+01  7.2080E-01  7.5070E-01  7.3466E-01
             9.9076E-01
 PARAMETER:  8.0799E-02 -4.3086E+00  6.0581E-02  5.6593E-01 -3.2478E-01  3.3249E-02  2.5371E+00 -2.2739E-01 -1.8676E-01 -2.0834E-01
             9.0717E-02
 GRADIENT:  -3.8451E-02 -1.4014E+00  1.1564E+00  3.8449E+00 -1.4065E+00  1.1473E-01 -2.9644E+00 -1.0878E-01  1.2647E+00  5.0607E-01
            -4.2711E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1684.41666012103        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1744
 NPARAMETR:  9.8091E-01  1.1942E-02  9.6664E-01  1.5942E+00  6.5575E-01  9.3528E-01  1.1605E+01  7.3767E-01  7.5177E-01  7.3139E-01
             9.9211E-01
 PARAMETER:  8.0726E-02 -4.3277E+00  6.6074E-02  5.6636E-01 -3.2198E-01  3.3091E-02  2.5515E+00 -2.0426E-01 -1.8533E-01 -2.1281E-01
             9.2083E-02
 GRADIENT:  -1.0163E-01 -1.2713E-01  3.6347E-01  4.2872E-01 -8.7859E-01 -1.8546E-02 -2.5229E-01  1.0965E-01  3.5445E-01  1.8875E-01
             2.4778E-01

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1684.41666398398        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1841
 NPARAMETR:  9.8093E-01  1.1951E-02  9.6699E-01  1.5943E+00  6.5595E-01  9.3531E-01  1.1596E+01  7.3771E-01  7.5167E-01  7.3139E-01
             9.9192E-01
 PARAMETER:  8.0741E-02 -4.3266E+00  6.6453E-02  5.6639E-01 -3.2166E-01  3.3116E-02  2.5509E+00 -2.0405E-01 -1.8545E-01 -2.1271E-01
             9.1990E-02
 GRADIENT:  -7.5874E-02  2.3051E+02  1.7112E-01 -1.7687E+03  1.8388E-02 -1.7731E-02  3.9052E+02  7.2039E-02  1.7243E-01  1.2439E-01
             1.8379E-01
 NUMSIGDIG:         3.5         3.3         3.0         3.3         5.2         3.5         3.3         2.3         3.7         2.5
                    2.2

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1841
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4526E-04  1.5067E-02 -1.7797E-02 -7.9458E-03 -1.8677E-02
 SE:             2.9821E-02  7.3098E-03  1.6130E-02  2.8800E-02  2.1323E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9344E-01  3.9282E-02  2.6988E-01  7.8263E-01  3.8106E-01

 ETASHRINKSD(%)  9.6022E-02  7.5511E+01  4.5962E+01  3.5159E+00  2.8566E+01
 ETASHRINKVR(%)  1.9195E-01  9.4003E+01  7.0799E+01  6.9081E+00  4.8972E+01
 EBVSHRINKSD(%)  4.5794E-01  8.0897E+01  4.6910E+01  3.4440E+00  2.6573E+01
 EBVSHRINKVR(%)  9.1379E-01  9.6351E+01  7.1815E+01  6.7694E+00  4.6085E+01
 RELATIVEINF(%)  9.9012E+01  2.9704E+00  2.9061E+00  6.6923E+01  5.5473E+00
 EPSSHRINKSD(%)  4.4321E+01
 EPSSHRINKVR(%)  6.8998E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1684.4166639839834     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -949.26583742024525     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.96
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1684.417       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.20E-02  9.67E-01  1.59E+00  6.56E-01  9.35E-01  1.16E+01  7.38E-01  7.52E-01  7.31E-01  9.92E-01
 


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
+        1.34E+03
 
 TH 2
+       -7.12E+02  9.35E+07
 
 TH 3
+       -1.28E+01  8.72E+03  1.79E+03
 
 TH 4
+        3.70E+01  1.60E+04 -6.72E+02  3.08E+05
 
 TH 5
+        7.62E+02 -2.22E+04 -1.23E+07  1.39E+03  5.64E+06
 
 TH 6
+        1.49E+01  9.03E+02  1.39E+02 -6.44E+01 -8.63E+02  2.16E+02
 
 TH 7
+       -1.19E+00 -6.25E+02  1.47E+01  2.70E+01 -3.74E+01  1.58E+00  2.86E+02
 
 TH 8
+       -1.21E+01  1.91E+03  2.30E+02 -1.42E+02 -7.90E+06  1.38E+01  3.29E+00  1.39E+02
 
 TH 9
+       -1.30E+02  1.25E+04  2.32E+03 -8.99E+02 -8.53E+06  1.88E+02  2.12E+01  6.10E+02  1.29E+07
 
 TH10
+       -3.54E+01  2.55E+03  6.07E+02 -2.19E+02 -7.64E+06  3.85E+01  4.28E+00  1.96E+02  1.27E+03  4.43E+02
 
 TH11
+        2.73E+01  4.97E+02  7.85E+01 -4.55E+01 -9.23E+02 -9.17E+00  8.52E-01  8.30E+01  1.76E+02  8.33E+01  2.49E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.743
Stop Time:
Sat Sep 25 10:05:06 CDT 2021
