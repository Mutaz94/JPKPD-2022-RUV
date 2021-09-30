Thu Sep 30 08:23:46 CDT 2021
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
$DATA ../../../../data/spa2/TD2/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2151.91590629309        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1737E+02 -4.7611E+01  2.0907E+01 -1.8268E+01  1.3295E+02  4.2688E+01 -2.4457E-02 -4.6304E+02 -1.2113E+02  1.8234E+00
            -1.0350E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2339.84933891295        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0073E+00  9.8856E-01  1.1345E+00  1.0118E+00  1.0283E+00  1.0403E+00  9.6666E-01  1.8905E+00  8.7823E-01  9.4532E-01
             9.6055E-01
 PARAMETER:  1.0723E-01  8.8497E-02  2.2621E-01  1.1171E-01  1.2793E-01  1.3953E-01  6.6092E-02  7.3684E-01 -2.9852E-02  4.3766E-02
             5.9747E-02
 GRADIENT:   5.1806E+02 -3.3735E+01 -1.9596E+01 -1.5234E+01  8.7919E+01  9.5228E+01 -8.5498E+00 -1.4107E+02 -4.1134E+01 -1.0020E+01
            -1.1117E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2359.56979657915        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      245
 NPARAMETR:  1.0051E+00  9.6697E-01  1.7225E+00  1.0861E+00  1.2039E+00  1.1157E+00  1.2907E+00  2.0787E+00  8.8453E-01  1.1158E+00
             1.0662E+00
 PARAMETER:  1.0511E-01  6.6413E-02  6.4376E-01  1.8261E-01  2.8561E-01  2.0948E-01  3.5516E-01  8.3173E-01 -2.2702E-02  2.0961E-01
             1.6413E-01
 GRADIENT:  -1.5420E+01 -4.0535E+01  7.3259E+00 -3.3967E+01  8.1462E+01  2.6774E+01  7.6913E+00 -1.7179E+02 -7.9453E+00 -1.3915E+00
             1.1622E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2401.17389045759        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      422
 NPARAMETR:  1.0167E+00  8.6867E-01  1.9445E+00  1.1767E+00  1.1498E+00  1.0457E+00  1.3629E+00  3.3145E+00  7.9904E-01  1.0244E+00
             1.0286E+00
 PARAMETER:  1.1658E-01 -4.0793E-02  7.6501E-01  2.6268E-01  2.3960E-01  1.4471E-01  4.0960E-01  1.2983E+00 -1.2434E-01  1.2410E-01
             1.2819E-01
 GRADIENT:   6.7167E+00 -1.7421E+01 -4.9764E+00  2.6167E+00  3.3272E+01  2.6092E+00 -2.5569E+00 -1.9192E+01 -8.1220E+00 -4.8727E+00
            -1.4671E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2404.48748799362        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      597
 NPARAMETR:  1.0121E+00  1.1168E+00  1.7843E+00  1.0306E+00  1.1951E+00  1.0366E+00  1.0550E+00  3.5916E+00  9.5836E-01  1.0843E+00
             1.0582E+00
 PARAMETER:  1.1205E-01  2.1047E-01  6.7903E-01  1.3016E-01  2.7820E-01  1.3593E-01  1.5353E-01  1.3786E+00  5.7464E-02  1.8091E-01
             1.5654E-01
 GRADIENT:  -4.7932E+00  4.6478E+00  2.4023E-01  1.3867E+01  1.5187E+00 -1.2985E+00 -1.1407E+00  6.3862E-01 -4.8351E-01 -1.9559E+00
             2.5403E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2404.66193775117        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      778
 NPARAMETR:  1.0128E+00  1.1395E+00  1.7632E+00  1.0069E+00  1.2019E+00  1.0409E+00  1.0475E+00  3.5810E+00  9.7027E-01  1.0987E+00
             1.0570E+00
 PARAMETER:  1.1271E-01  2.3059E-01  6.6712E-01  1.0686E-01  2.8391E-01  1.4012E-01  1.4636E-01  1.3757E+00  6.9823E-02  1.9411E-01
             1.5543E-01
 GRADIENT:  -3.3756E+00 -2.0575E+00  1.1998E+00  1.9118E-01 -7.4619E-01  3.6963E-01  1.3949E-01 -4.5100E-01  5.5922E-02 -1.5594E+00
             3.5992E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2404.68552866474        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      958             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0137E+00  1.1402E+00  1.7293E+00  1.0074E+00  1.2045E+00  1.0414E+00  1.0470E+00  3.5492E+00  9.7000E-01  1.1055E+00
             1.0573E+00
 PARAMETER:  1.1357E-01  2.3117E-01  6.4769E-01  1.0733E-01  2.8604E-01  1.4058E-01  1.4589E-01  1.3667E+00  6.9543E-02  2.0028E-01
             1.5569E-01
 GRADIENT:   4.7121E+02  1.5134E+02  1.5767E+01  8.0499E+01  5.4868E+01  8.3380E+01  1.6899E+01  4.4081E+01  6.8412E+00  3.9942E+00
             2.9442E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2404.71154004425        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1139             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0155E+00  1.1445E+00  1.7329E+00  1.0021E+00  1.2045E+00  1.0412E+00  1.0429E+00  3.5446E+00  9.7405E-01  1.1088E+00
             1.0557E+00
 PARAMETER:  1.1542E-01  2.3496E-01  6.4979E-01  1.0207E-01  2.8607E-01  1.4039E-01  1.4198E-01  1.3654E+00  7.3706E-02  2.0325E-01
             1.5425E-01
 GRADIENT:   4.8428E+02  1.5394E+02  1.6738E+01  7.2354E+01  5.2857E+01  8.3022E+01  1.6493E+01  4.3712E+01  7.0168E+00  4.3577E+00
             1.1871E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2404.72084919233        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1321
 NPARAMETR:  1.0148E+00  1.1504E+00  1.7389E+00  1.0012E+00  1.2039E+00  1.0404E+00  1.0403E+00  3.5380E+00  9.7513E-01  1.1122E+00
             1.0565E+00
 PARAMETER:  1.1471E-01  2.4014E-01  6.5325E-01  1.0122E-01  2.8561E-01  1.3960E-01  1.3955E-01  1.3636E+00  7.4814E-02  2.0635E-01
             1.5492E-01
 GRADIENT:   7.0378E-01 -3.3966E-01  1.0130E+00  2.5211E+00 -8.8965E-01  1.5281E-01  1.4393E-01 -2.9931E+00 -2.8803E-02 -1.2503E-01
            -2.5824E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2404.73856425917        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1500
 NPARAMETR:  1.0114E+00  1.1779E+00  1.7386E+00  9.8140E-01  1.2178E+00  1.0373E+00  1.0169E+00  3.5399E+00  9.8976E-01  1.1259E+00
             1.0594E+00
 PARAMETER:  1.1136E-01  2.6370E-01  6.5307E-01  8.1224E-02  2.9706E-01  1.3665E-01  1.1674E-01  1.3641E+00  8.9712E-02  2.1856E-01
             1.5766E-01
 GRADIENT:  -6.4689E+00 -1.8801E+00  2.0827E+00 -9.5674E-01 -6.0378E-01 -1.0899E+00 -1.7387E-01 -4.0834E+00 -8.8569E-02  1.7210E-02
            -4.5770E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2404.76969051459        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1682
 NPARAMETR:  1.0158E+00  1.1813E+00  1.6988E+00  9.8047E-01  1.2149E+00  1.0412E+00  1.0178E+00  3.5040E+00  9.9077E-01  1.1269E+00
             1.0590E+00
 PARAMETER:  1.1564E-01  2.6660E-01  6.2995E-01  8.0272E-02  2.9463E-01  1.4033E-01  1.1767E-01  1.3539E+00  9.0723E-02  2.1948E-01
             1.5733E-01
 GRADIENT:   2.4217E+00 -3.1526E-01  6.6970E-01  1.5972E+00  4.1094E-02  4.0944E-01  2.2891E-01 -5.4564E+00 -1.8699E-01  1.2063E-01
             1.9833E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -2404.76969051459        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1709
 NPARAMETR:  1.0157E+00  1.1808E+00  1.7004E+00  9.7990E-01  1.2148E+00  1.0412E+00  1.0166E+00  3.4970E+00  9.9131E-01  1.1267E+00
             1.0589E+00
 PARAMETER:  1.1564E-01  2.6660E-01  6.2995E-01  8.0272E-02  2.9463E-01  1.4033E-01  1.1767E-01  1.3539E+00  9.0723E-02  2.1948E-01
             1.5733E-01
 GRADIENT:   1.7115E-01  5.2235E+03 -2.2070E+03  1.3294E+00  2.2589E-01 -7.7588E-02  2.0789E-01  2.7221E+03 -1.5388E-01  1.0490E-01
             1.9959E-01
 NUMSIGDIG:         3.1         2.3         2.3         1.7         3.1         2.9         1.5         2.3         1.7         2.6
                    2.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1709
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0746E-04 -2.9231E-02 -4.3520E-02  1.5196E-02 -5.8912E-02
 SE:             2.9890E-02  1.9818E-02  2.3060E-02  2.4451E-02  2.1281E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8379E-01  1.4022E-01  5.9126E-02  5.3427E-01  5.6352E-03

 ETASHRINKSD(%)  1.0000E-10  3.3606E+01  2.2747E+01  1.8087E+01  2.8705E+01
 ETASHRINKVR(%)  1.0000E-10  5.5919E+01  4.0319E+01  3.2902E+01  4.9171E+01
 EBVSHRINKSD(%)  3.4328E-01  3.3427E+01  2.3070E+01  2.0543E+01  2.6465E+01
 EBVSHRINKVR(%)  6.8538E-01  5.5681E+01  4.0818E+01  3.6866E+01  4.5927E+01
 RELATIVEINF(%)  9.9295E+01  9.7515E+00  3.6361E+01  1.4892E+01  3.0827E+01
 EPSSHRINKSD(%)  3.1998E+01
 EPSSHRINKVR(%)  5.3757E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2404.7696905145917     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1302.0434506689846     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.17
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.99
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2404.770       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.18E+00  1.70E+00  9.80E-01  1.21E+00  1.04E+00  1.02E+00  3.50E+00  9.91E-01  1.13E+00  1.06E+00
 


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
+        2.53E+06
 
 TH 2
+       -4.73E+00  5.99E+05
 
 TH 3
+       -1.88E+01  2.99E+02  3.03E+04
 
 TH 4
+        2.07E+02  1.93E+06  3.86E+02  7.19E+02
 
 TH 5
+        8.29E+05 -5.28E+05 -6.19E+01 -9.93E+05  7.37E+05
 
 TH 6
+        5.60E-01  3.14E+01 -9.49E+00 -1.24E+00 -3.77E-01  1.82E+02
 
 TH 7
+        7.62E-01 -9.25E+05 -1.01E+02 -1.47E+01 -8.13E+00 -2.05E-01  5.13E+01
 
 TH 8
+       -6.08E+04  6.17E+04  4.25E+02  7.27E+04 -5.43E+04  5.47E+00 -5.96E+04  4.02E+03
 
 TH 9
+        1.01E+00 -1.50E+01 -1.04E+02  1.65E+01 -9.83E+05 -1.70E-01  3.41E+01  7.20E+04  8.64E+01
 
 TH10
+       -8.36E+01 -1.46E+01 -1.38E+01  1.44E+06  3.94E+05 -1.28E-01 -5.44E+00 -2.89E+04  5.40E+00  5.70E+05
 
 TH11
+       -6.40E+00 -2.71E+01 -2.67E+01  4.27E+06  5.84E+05  2.37E+00  8.28E+00 -4.29E+04  2.79E+00  8.46E+05  4.54E+02
 
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
 #CPUT: Total CPU Time in Seconds,       49.239
Stop Time:
Thu Sep 30 08:24:37 CDT 2021
