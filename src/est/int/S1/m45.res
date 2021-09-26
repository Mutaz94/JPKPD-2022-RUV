Fri Sep 24 23:21:00 CDT 2021
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
$DATA ../../../../data/int/S1/dat45.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3360.81736754964        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0870E+02 -3.1300E+01  1.1001E+02  1.2404E+00  1.2089E+02 -6.4170E+01 -2.8699E+01 -4.9233E+02 -1.2770E+02 -2.1193E+01
            -2.9488E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3585.44993853426        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7178E-01  9.7834E-01  9.7276E-01  9.9273E-01  9.5755E-01  1.0984E+00  1.0327E+00  2.0353E+00  1.0149E+00  9.8758E-01
             1.0213E+00
 PARAMETER:  7.1370E-02  7.8102E-02  7.2379E-02  9.2705E-02  5.6626E-02  1.9390E-01  1.3221E-01  8.1063E-01  1.1475E-01  8.7502E-02
             1.2111E-01
 GRADIENT:   4.3277E+01 -2.3352E+01  5.5092E+00 -1.9099E+01  4.9631E+01 -1.0439E+01 -2.3838E+01 -1.1954E+02 -6.6060E+00 -1.3068E+01
            -1.4268E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3587.67262734273        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      179
 NPARAMETR:  1.0065E+00  9.9748E-01  1.0348E+00  9.7045E-01  9.7852E-01  1.0978E+00  1.0654E+00  2.0756E+00  1.0441E+00  1.0226E+00
             1.0329E+00
 PARAMETER:  1.0650E-01  9.7475E-02  1.3425E-01  7.0009E-02  7.8290E-02  1.9329E-01  1.6338E-01  8.3026E-01  1.4316E-01  1.2235E-01
             1.3241E-01
 GRADIENT:   1.0960E+02 -2.7357E+01  1.7390E+01 -3.8920E+01  3.2221E+01 -1.4355E+01 -1.8083E+01 -1.1318E+02 -1.0537E-01 -4.5012E+00
            -1.2176E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3588.66535556736        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  1.0063E+00  9.9731E-01  1.0348E+00  9.7070E-01  9.7843E-01  1.1831E+00  1.0664E+00  2.0779E+00  1.0494E+00  1.0225E+00
             1.0330E+00
 PARAMETER:  1.0624E-01  9.7302E-02  1.3423E-01  7.0259E-02  7.8198E-02  2.6811E-01  1.6425E-01  8.3138E-01  1.4817E-01  1.2228E-01
             1.3244E-01
 GRADIENT:   1.0021E+02 -2.7502E+01  1.7352E+01 -3.8093E+01  3.2138E+01  2.1790E+01 -1.7930E+01 -1.1251E+02  1.1447E+00 -4.3750E+00
            -1.2112E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3588.68943580157        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0063E+00  9.9731E-01  1.0348E+00  9.7070E-01  9.7843E-01  1.1703E+00  1.0664E+00  2.0779E+00  1.0495E+00  1.0225E+00
             1.0330E+00
 PARAMETER:  1.0624E-01  9.7302E-02  1.3423E-01  7.0259E-02  7.8198E-02  2.5726E-01  1.6425E-01  8.3137E-01  1.4827E-01  1.2228E-01
             1.3244E-01
 GRADIENT:   5.6843E+01 -3.5256E+01  1.6461E+01 -4.7173E+01  2.5975E+01  1.2076E-01 -1.9220E+01 -1.1631E+02 -3.8662E-02 -8.0806E+00
            -1.2140E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3595.59603932063        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.0059E+00  9.9753E-01  1.0346E+00  9.7099E-01  9.7827E-01  1.1357E+00  1.0667E+00  2.3339E+00  1.2488E+00  1.0226E+00
             1.0344E+00
 PARAMETER:  1.0584E-01  9.7526E-02  1.3404E-01  7.0558E-02  7.8033E-02  2.2726E-01  1.6458E-01  9.4754E-01  3.2216E-01  1.2236E-01
             1.3380E-01
 GRADIENT:   5.9763E+01 -3.6447E+01  1.1603E+01 -3.0170E+01  2.1555E+01 -1.2131E+01 -1.9421E+01 -6.7608E+01  4.1922E+01 -7.1003E+00
            -1.0525E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3598.89434595514        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:      885             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0057E+00  9.9860E-01  1.0346E+00  9.7187E-01  9.7822E-01  1.1384E+00  1.0668E+00  2.4208E+00  1.2225E+00  1.0237E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  9.8601E-02  1.3399E-01  7.1462E-02  7.7981E-02  2.2961E-01  1.6469E-01  9.8412E-01  3.0092E-01  1.2347E-01
             1.3422E-01
 GRADIENT:   1.0345E+02 -2.5513E+01  1.0307E+01 -2.1818E+01  2.5635E+01  3.7229E+00 -1.7999E+01 -5.2592E+01  4.2236E+01 -3.4903E+00
            -1.0223E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3599.29182018037        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      996
 NPARAMETR:  1.0057E+00  9.9922E-01  1.0346E+00  9.7280E-01  9.7822E-01  1.1380E+00  1.0668E+00  2.4208E+00  1.2117E+00  1.0249E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  9.9225E-02  1.3399E-01  7.2419E-02  7.7981E-02  2.2928E-01  1.6469E-01  9.8412E-01  2.9206E-01  1.2455E-01
             1.3422E-01
 GRADIENT:   5.9273E+01 -3.2039E+01  9.3558E+00 -3.0190E+01  1.9259E+01 -1.1291E+01 -1.9168E+01 -5.7227E+01  3.7948E+01 -7.3568E+00
            -1.0263E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3645.55622317089        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1174
 NPARAMETR:  1.0057E+00  1.0022E+00  1.0346E+00  9.7760E-01  9.7822E-01  1.1462E+00  1.0668E+00  2.4208E+00  1.1705E+00  1.0274E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  1.0219E-01  1.3398E-01  7.7341E-02  7.7979E-02  2.3643E-01  1.6468E-01  9.8410E-01  2.5741E-01  1.2701E-01
             1.3422E-01
 GRADIENT:   1.0474E+02 -4.9895E+01  6.2087E+00 -1.4097E+01 -1.7818E+01  9.0736E+00  1.3422E+00 -3.9265E+01  3.3569E+01  3.2875E+01
            -9.0207E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3645.73426979766        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:     1327
 NPARAMETR:  1.0057E+00  1.0022E+00  1.0346E+00  9.7760E-01  9.7824E-01  1.1693E+00  1.0669E+00  2.4214E+00  1.1671E+00  1.0274E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  1.0217E-01  1.3402E-01  7.7341E-02  7.8004E-02  2.5645E-01  1.6472E-01  9.8435E-01  2.5456E-01  1.2704E-01
             1.3419E-01
 GRADIENT:   5.6147E+01 -5.8817E+01  5.2970E+00 -2.3703E+01 -2.2674E+01 -1.9747E-01 -7.2465E-01 -4.2226E+01  3.0582E+01  3.2340E+01
            -9.0492E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3646.21947964096        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1471
 NPARAMETR:  1.0057E+00  1.0022E+00  1.0346E+00  9.7760E-01  9.7824E-01  1.1690E+00  1.0669E+00  2.4377E+00  1.1591E+00  1.0274E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  1.0217E-01  1.3402E-01  7.7341E-02  7.8004E-02  2.5617E-01  1.6472E-01  9.9104E-01  2.4763E-01  1.2703E-01
             1.3419E-01
 GRADIENT:   1.0252E+02 -4.9469E+01  5.7999E+00 -1.5014E+01 -1.8121E+01  1.8517E+01  1.2861E+00 -3.7607E+01  3.1701E+01  3.2710E+01
            -8.9954E+01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -3646.21947964096        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1546
 NPARAMETR:  1.0057E+00  1.0022E+00  1.0346E+00  9.7760E-01  9.7824E-01  1.1693E+00  1.0668E+00  2.4375E+00  1.1591E+00  1.0274E+00
             1.0348E+00
 PARAMETER:  1.0571E-01  1.0217E-01  1.3402E-01  7.7341E-02  7.8004E-02  2.5617E-01  1.6472E-01  9.9104E-01  2.4763E-01  1.2703E-01
             1.3419E-01
 GRADIENT:   5.6570E+01 -2.9658E+05  2.2607E+05 -4.2291E+01  7.4792E+00 -3.0334E-01  3.6783E+05  6.1034E+04  1.2237E+05  2.3851E+05
            -2.2587E+05
 NUMSIGDIG:         7.6         3.3         3.3         7.8         8.5         2.2         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1546
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.9138E-02 -1.2692E-02 -4.1004E-02  3.2515E-02 -4.6812E-02
 SE:             2.9797E-02  2.2186E-02  2.8593E-02  2.4658E-02  2.0573E-02
 N:                     100         100         100         100         100

 P VAL.:         3.2814E-01  5.6729E-01  1.5155E-01  1.8729E-01  2.2885E-02

 ETASHRINKSD(%)  1.7671E-01  2.5673E+01  4.2102E+00  1.7394E+01  3.1076E+01
 ETASHRINKVR(%)  3.5311E-01  4.4755E+01  8.2430E+00  3.1762E+01  5.2495E+01
 EBVSHRINKSD(%)  2.1545E-01  2.4773E+01  1.6259E+01  8.4825E+00  2.3419E+01
 EBVSHRINKVR(%)  4.3044E-01  4.3410E+01  2.9874E+01  1.6246E+01  4.1354E+01
 RELATIVEINF(%)  9.9569E+01  3.1726E+01  6.1295E+01  6.1721E+01  3.2311E+01
 EPSSHRINKSD(%)  2.0279E+01
 EPSSHRINKVR(%)  3.6446E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3646.2194796409558     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1992.1301198725450     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    48.55
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3646.219       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E+00  1.03E+00  9.78E-01  9.78E-01  1.17E+00  1.07E+00  2.44E+00  1.16E+00  1.03E+00  1.03E+00
 


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
+        2.68E+09
 
 TH 2
+        3.61E+00  1.44E+09
 
 TH 3
+        1.03E+09  3.39E+04  1.58E+09
 
 TH 4
+        1.46E+09  2.31E+02 -1.58E+01  3.17E+09
 
 TH 5
+        1.46E+09 -1.51E+09  1.12E+09  1.58E+09  3.17E+09
 
 TH 6
+       -4.88E+01 -1.15E+03  8.49E+02 -1.01E+00 -2.20E+00  1.45E+02
 
 TH 7
+        8.11E+08  2.68E+04  6.22E+08  1.04E+01  8.81E+08  6.68E+02  9.81E+08
 
 TH 8
+       -8.26E-02  1.31E+05 -9.66E+04  1.01E+00  6.40E+07  4.86E+01 -7.62E+04  2.59E+06
 
 TH 9
+        3.10E+00  1.15E+04 -8.51E+03  4.46E+01 -2.32E+00  4.09E+02 -6.71E+03 -4.80E+02  1.84E+08
 
 TH10
+        1.09E+09  3.60E+04  8.37E+08  1.43E+01  1.19E+09  9.01E+02 -5.85E+04 -1.03E+05 -9.03E+03  8.89E+08
 
 TH11
+       -1.03E+09  1.07E+09 -7.87E+08 -1.12E+09 -1.12E+09 -8.45E+02 -6.21E+08 -4.51E+07  8.50E+03 -8.36E+08  1.57E+09
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       63.767
Stop Time:
Fri Sep 24 23:22:06 CDT 2021
