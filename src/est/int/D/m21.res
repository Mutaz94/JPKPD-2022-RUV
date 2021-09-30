Wed Sep 29 08:13:40 CDT 2021
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
$DATA ../../../../data/int/D/dat21.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1062.34952561153        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9312E+01  1.1053E+02 -4.5663E+01 -1.9764E+02  2.6998E+02 -8.6980E+02 -2.8702E+02 -7.8752E+01 -8.5000E+02 -4.3253E+02
            -7.5903E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2374.61031019319        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      107
 NPARAMETR:  1.1578E+00  8.8296E-01  1.0476E+00  1.3284E+00  7.3999E-01  2.7138E+00  1.5116E+00  1.0256E+00  2.8458E+00  2.3353E+00
             4.0029E+00
 PARAMETER:  2.4649E-01 -2.4476E-02  1.4648E-01  3.8394E-01 -2.0112E-01  1.0983E+00  5.1316E-01  1.2532E-01  1.1458E+00  9.4815E-01
             1.4870E+00
 GRADIENT:  -2.3216E+01 -5.7256E+00  4.6069E+00  1.3866E+01 -3.5083E+01  4.9525E+01  1.9207E+01  1.1350E+01  5.1815E+01  1.0718E+01
             5.6219E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2507.27447732671        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      289
 NPARAMETR:  1.2503E+00  2.8390E+00  4.3682E+00  6.4312E-01  2.4314E+00  3.1475E+00  3.9282E+00  6.5739E-01  5.3979E+00  1.5909E+00
             3.4837E+00
 PARAMETER:  3.2340E-01  1.1434E+00  1.5744E+00 -3.4142E-01  9.8847E-01  1.2466E+00  1.4682E+00 -3.1948E-01  1.7860E+00  5.6429E-01
             1.3481E+00
 GRADIENT:   3.0892E+00  1.7636E+01 -1.9626E+01  2.8044E+01  1.5319E+01  9.3797E+01  7.4811E+01  4.7832E-01  5.3327E+01  3.4689E+01
             3.9028E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2574.53445473571        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      466
 NPARAMETR:  1.4660E+00  2.8020E+00  1.5926E+01  2.4136E-01  2.4699E+00  3.0280E+00  2.5491E+00  3.2479E+00  9.2300E+00  1.3472E+00
             2.8391E+00
 PARAMETER:  4.8252E-01  1.1303E+00  2.8680E+00 -1.3214E+00  1.0042E+00  1.2079E+00  1.0357E+00  1.2780E+00  2.3225E+00  3.9801E-01
             1.1435E+00
 GRADIENT:   4.8679E+01 -2.4775E+01 -1.5534E+01  1.1513E+01  1.1389E+01  7.3075E+01  4.1344E+01 -2.2321E-01  4.1927E+01 -7.0973E-01
            -1.5335E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2613.14750332388        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      643
 NPARAMETR:  1.1710E+00  2.7583E+00  8.5786E+01  1.1120E-01  2.4644E+00  2.3682E+00  1.8823E+00  1.6108E+00  8.5428E+00  1.2966E+00
             2.8594E+00
 PARAMETER:  2.5786E-01  1.1146E+00  4.5519E+00 -2.0964E+00  1.0020E+00  9.6211E-01  7.3250E-01  5.7675E-01  2.2451E+00  3.5974E-01
             1.1506E+00
 GRADIENT:  -1.7716E+01  2.7189E+01 -2.2843E+00 -8.9244E+00  2.5611E+01  3.2554E+00  2.4140E+01  5.1779E-01 -2.4610E+01 -1.2720E+01
             4.8535E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2614.97254673132        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      819
 NPARAMETR:  1.2273E+00  2.7091E+00  1.1547E+02  8.5318E-02  2.4398E+00  2.3366E+00  1.7894E+00  1.3557E+00  9.0744E+00  1.3522E+00
             2.8470E+00
 PARAMETER:  3.0478E-01  1.0966E+00  4.8490E+00 -2.3614E+00  9.9190E-01  9.4868E-01  6.8189E-01  4.0434E-01  2.3055E+00  4.0171E-01
             1.1463E+00
 GRADIENT:   1.4411E+00  2.7693E+01 -7.7159E-01 -9.0639E+00  1.0691E+01 -2.2564E+00  2.1872E+01  5.9077E-01 -2.8070E+01 -2.7129E+00
            -3.6167E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2626.21817882276        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  1.2493E+00  2.6503E+00  1.5471E+02  6.5735E-02  2.4209E+00  2.3492E+00  1.7252E+00  1.0945E+00  9.3009E+00  1.3860E+00
             2.8451E+00
 PARAMETER:  3.2256E-01  1.0747E+00  5.1415E+00 -2.6221E+00  9.8413E-01  9.5406E-01  6.4535E-01  1.9025E-01  2.3301E+00  4.2640E-01
             1.1456E+00
 GRADIENT:   8.9041E+00  1.4481E+01  4.9786E-01 -1.0872E+01  4.9525E-01 -2.7478E-01  1.9583E+01  3.4544E-01 -4.8880E+01  3.2285E+00
            -3.0397E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2626.51948381226        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  1.2247E+00  2.6516E+00  1.3233E+02  6.5073E-02  2.4182E+00  2.3493E+00  1.7258E+00  9.1930E-01  9.3661E+00  1.3665E+00
             2.8474E+00
 PARAMETER:  3.0267E-01  1.0752E+00  4.9853E+00 -2.6322E+00  9.8304E-01  9.5414E-01  6.4567E-01  1.5855E-02  2.3371E+00  4.1224E-01
             1.1464E+00
 GRADIENT:   6.3886E-01  1.4056E+01 -5.9519E-03 -1.0248E+01  1.8246E+00  1.7884E-01  1.8769E+01  3.0331E-01 -4.8522E+01  3.4735E-01
            -1.5316E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2629.41999347781        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1394             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2203E+00  2.6652E+00  1.3382E+02  6.2562E-02  2.4055E+00  2.3394E+00  1.7310E+00  8.3574E-01  9.4701E+00  1.3612E+00
             2.8351E+00
 PARAMETER:  2.9914E-01  1.0803E+00  4.9965E+00 -2.6716E+00  9.7777E-01  9.4991E-01  6.4869E-01 -7.9438E-02  2.3481E+00  4.0839E-01
             1.1421E+00
 GRADIENT:   1.8776E+02  4.5921E+02  4.9991E-01  4.6494E+01  2.6885E+01  2.1301E+02  7.5134E+01  2.5497E-01  1.0597E+03  1.8122E+00
             9.0872E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2629.45993477037        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1576
 NPARAMETR:  1.2204E+00  2.6648E+00  1.3378E+02  6.2692E-02  2.4059E+00  2.3394E+00  1.7307E+00  8.1688E-01  9.4715E+00  1.3613E+00
             2.8369E+00
 PARAMETER:  2.9914E-01  1.0801E+00  4.9962E+00 -2.6695E+00  9.7793E-01  9.4991E-01  6.4854E-01 -1.0227E-01  2.3483E+00  4.0841E-01
             1.1427E+00
 GRADIENT:  -7.1312E-01  5.2156E+00  1.1301E-01 -8.7664E+00 -2.6897E+00 -1.5178E+00  1.0676E+01  2.4412E-01 -4.9880E+01 -5.2525E-01
            -1.0112E+01

0ITERATION NO.:   49    OBJECTIVE VALUE:  -2629.46522865838        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:     1720
 NPARAMETR:  1.2204E+00  2.6647E+00  1.3377E+02  6.2709E-02  2.4060E+00  2.3394E+00  1.7307E+00  8.1259E-01  9.4718E+00  1.3613E+00
             2.8371E+00
 PARAMETER:  2.9914E-01  1.0801E+00  4.9961E+00 -2.6693E+00  9.7795E-01  9.4991E-01  6.4853E-01 -1.0647E-01  2.3483E+00  4.0842E-01
             1.1428E+00
 GRADIENT:  -2.6441E+00  2.8860E+03  4.8582E+00 -1.1855E+03 -2.3179E+00  4.0616E+01  4.7710E+03  2.4102E-01 -1.2153E+03 -2.9962E-01
             3.1986E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1720
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7833E-03 -1.7772E-02 -1.9719E-04  3.8783E-02 -1.5170E-02
 SE:             3.0544E-02  2.7140E-02  4.0279E-04  1.5374E-02  2.5661E-02
 N:                     100         100         100         100         100

 P VAL.:         9.0142E-01  5.1258E-01  6.2444E-01  1.1648E-02  5.5440E-01

 ETASHRINKSD(%)  1.0000E-10  9.0791E+00  9.8651E+01  4.8495E+01  1.4032E+01
 ETASHRINKVR(%)  1.0000E-10  1.7334E+01  9.9982E+01  7.3473E+01  2.6095E+01
 EBVSHRINKSD(%)  3.4207E-01  6.4942E+00  9.6284E+01  6.6908E+01  1.3349E+01
 EBVSHRINKVR(%)  6.8296E-01  1.2567E+01  9.9862E+01  8.9049E+01  2.4916E+01
 RELATIVEINF(%)  9.9310E+01  6.2265E+01  1.3454E-01  7.8074E+00  7.3039E+01
 EPSSHRINKSD(%)  1.5349E+01
 EPSSHRINKVR(%)  2.8341E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2629.4652286583828     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -975.37586888997203     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    64.53
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2629.465       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.22E+00  2.66E+00  1.34E+02  6.27E-02  2.41E+00  2.34E+00  1.73E+00  8.13E-01  9.47E+00  1.36E+00  2.84E+00
 


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
+        1.35E+02
 
 TH 2
+        1.74E+00  1.03E+04
 
 TH 3
+       -3.82E-04 -8.60E-03  2.80E-03
 
 TH 4
+       -3.66E+01 -4.72E+02  1.45E-01  3.36E+06
 
 TH 5
+       -6.89E-01 -2.74E+01 -1.65E-02  4.28E+02  2.95E+02
 
 TH 6
+        1.48E+03  4.38E+00 -2.78E-04 -7.54E+01 -2.30E+02  2.80E+02
 
 TH 7
+        5.49E+00  2.54E+04 -2.18E-02  8.78E+02 -6.46E+01  1.08E+01  6.61E+04
 
 TH 8
+       -1.08E-01 -2.58E+00 -1.14E-03  4.42E+01 -2.89E-02 -6.23E+03 -6.60E+00  8.96E-01
 
 TH 9
+       -2.75E-01 -1.11E+02  1.45E-03  3.48E+03  3.13E+00 -5.83E-01 -2.72E+02  3.38E-01  2.03E+02
 
 TH10
+       -6.29E-02 -1.71E+01 -8.68E-04  2.93E+02 -9.25E+02  9.71E+02 -4.37E+01  5.81E-03  2.15E+00  3.94E+03
 
 TH11
+       -3.05E+00 -1.23E+01 -2.86E-03  1.60E+02 -1.57E+02  1.67E+02 -2.14E+01 -3.38E-02  1.10E+00  6.71E+02  2.59E+02
 
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
 #CPUT: Total CPU Time in Seconds,       84.276
Stop Time:
Wed Sep 29 08:15:06 CDT 2021
