Wed Sep 29 21:54:02 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat19.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1398.32710315537        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7266E+02  2.1089E+01  4.2845E+01  3.0848E+01  5.9955E+01  6.0227E+01 -3.6520E+00  5.5578E-01  4.3284E+01 -4.8461E+01
            -1.4729E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1832.27912913891        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0288E+00  1.0918E+00  1.0165E+00  1.0245E+00  1.0243E+00  8.6024E-01  9.2366E-01  8.7732E-01  6.7286E-01  1.0406E+00
             1.9930E+00
 PARAMETER:  1.2838E-01  1.8782E-01  1.1640E-01  1.2416E-01  1.2400E-01 -5.0543E-02  2.0592E-02 -3.0881E-02 -2.9622E-01  1.3976E-01
             7.8965E-01
 GRADIENT:   1.0557E+02  6.4723E+01 -6.2851E-01  9.2790E+01 -6.9784E-01 -2.5724E+01 -1.8139E+01  4.8381E+00 -1.5487E+01 -5.9392E+00
            -7.0056E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1837.01627413075        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0186E+00  9.9266E-01  8.4600E-01  1.0652E+00  9.0130E-01  8.6328E-01  9.9043E-01  2.6117E-01  7.3547E-01  9.0660E-01
             2.0087E+00
 PARAMETER:  1.1842E-01  9.2629E-02 -6.7237E-02  1.6320E-01 -3.9130E-03 -4.7019E-02  9.0385E-02 -1.2426E+00 -2.0724E-01  1.9425E-03
             7.9750E-01
 GRADIENT:   6.1655E+01  3.1700E+01 -7.1019E+00  8.7698E+01  2.4581E+01 -2.5236E+01 -1.0545E+01  4.9923E-01 -2.4380E+00 -1.2458E+01
            -5.2937E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1841.09452263119        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0166E+00  9.8265E-01  7.6056E-01  1.0332E+00  8.4367E-01  9.2138E-01  1.0623E+00  2.2933E-01  7.1963E-01  9.0608E-01
             2.0592E+00
 PARAMETER:  1.1643E-01  8.2497E-02 -1.7370E-01  1.3269E-01 -6.9992E-02  1.8115E-02  1.6039E-01 -1.3726E+00 -2.2902E-01  1.3692E-03
             8.2230E-01
 GRADIENT:  -6.3199E+01 -6.7866E+00 -2.2967E+00 -1.4392E+01  8.8312E+00 -7.8963E+00 -6.0026E+00  5.3149E-01 -2.1661E+00 -1.7873E+00
            -2.5396E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1845.17315019440        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      479
 NPARAMETR:  1.0383E+00  6.7823E-01  5.2883E-01  1.1747E+00  5.6271E-01  9.2911E-01  1.4585E+00  6.0129E-02  6.8963E-01  6.3551E-01
             2.1049E+00
 PARAMETER:  1.3755E-01 -2.8827E-01 -5.3708E-01  2.6103E-01 -4.7499E-01  2.6474E-02  4.7742E-01 -2.7113E+00 -2.7160E-01 -3.5333E-01
             8.4425E-01
 GRADIENT:  -1.3403E+01  1.1208E+01 -7.7239E+00  3.7790E+01  3.4217E+00 -3.8669E+00 -6.3172E-01  5.7074E-02  6.9310E-01 -6.3952E-01
             6.3640E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1845.68424257766        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  1.0445E+00  6.1528E-01  5.4647E-01  1.1929E+00  5.5533E-01  9.3422E-01  1.5666E+00  4.2025E-02  6.8608E-01  6.5888E-01
             2.0831E+00
 PARAMETER:  1.4353E-01 -3.8567E-01 -5.0428E-01  2.7637E-01 -4.8820E-01  3.1962E-02  5.4893E-01 -3.0695E+00 -2.7676E-01 -3.1722E-01
             8.3388E-01
 GRADIENT:   6.3250E+00 -1.1711E-01  1.2401E+00 -3.6918E+00 -2.8843E+00 -8.2047E-01 -9.4363E-01  2.8968E-02  2.5560E+00  3.3308E-01
            -6.9111E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1846.14584739209        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0357E+00  4.7867E-01  7.6991E-01  1.3155E+00  6.5680E-01  9.3220E-01  1.9813E+00  1.1896E-02  6.2816E-01  8.3606E-01
             2.0915E+00
 PARAMETER:  1.3506E-01 -6.3675E-01 -1.6148E-01  3.7418E-01 -3.2037E-01  2.9789E-02  7.8377E-01 -4.3315E+00 -3.6497E-01 -7.9054E-02
             8.3789E-01
 GRADIENT:  -1.9577E-01  2.9603E+00  3.2331E-01  7.5783E+00 -2.0939E+00  9.7585E-01  3.5393E-01  1.6954E-03 -8.4187E-01  2.6458E-01
            -2.8296E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1846.20816150921        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1012             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0334E+00  3.9500E-01  8.4392E-01  1.3679E+00  6.7917E-01  9.2747E-01  2.2351E+00  1.0000E-02  6.2029E-01  8.8707E-01
             2.0980E+00
 PARAMETER:  1.3284E-01 -8.2886E-01 -6.9696E-02  4.1330E-01 -2.8689E-01  2.4708E-02  9.0426E-01 -5.3711E+00 -3.7758E-01 -1.9833E-02
             8.4099E-01
 GRADIENT:   1.2709E+02  1.1950E+01  4.1581E-01  1.4052E+02  9.7774E+00  7.9610E+00  6.5643E+00  0.0000E+00  6.3889E+00  4.3027E-01
             8.2877E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1846.20816150921        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1076
 NPARAMETR:  1.0333E+00  3.9596E-01  8.4477E-01  1.3681E+00  6.7859E-01  9.2739E-01  2.2398E+00  1.0000E-02  6.1983E-01  8.8665E-01
             2.0982E+00
 PARAMETER:  1.3284E-01 -8.2886E-01 -6.9696E-02  4.1330E-01 -2.8689E-01  2.4708E-02  9.0426E-01 -5.3711E+00 -3.7758E-01 -1.9833E-02
             8.4099E-01
 GRADIENT:   1.4855E-01 -3.9291E-01 -4.3674E-01 -4.1366E-01  1.3201E+00  3.8005E-02 -2.2415E-01  0.0000E+00  1.4017E-01  3.0527E-02
            -9.3738E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1076
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6124E-03  2.8677E-02 -1.4145E-04 -2.7973E-02 -8.0712E-03
 SE:             2.9437E-02  1.8151E-02  1.5996E-04  2.3019E-02  2.0626E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5632E-01  1.1414E-01  3.7655E-01  2.2429E-01  6.9556E-01

 ETASHRINKSD(%)  1.3837E+00  3.9191E+01  9.9464E+01  2.2883E+01  3.0901E+01
 ETASHRINKVR(%)  2.7483E+00  6.3023E+01  9.9997E+01  4.0529E+01  5.2253E+01
 EBVSHRINKSD(%)  1.5566E+00  4.4386E+01  9.9413E+01  2.0412E+01  2.8666E+01
 EBVSHRINKVR(%)  3.0890E+00  6.9071E+01  9.9997E+01  3.6657E+01  4.9114E+01
 RELATIVEINF(%)  9.5697E+01  4.3152E+00  2.3323E-04  8.8611E+00  3.7897E+00
 EPSSHRINKSD(%)  2.7135E+01
 EPSSHRINKVR(%)  4.6906E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1846.2081615092102     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.26962830453749     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.98
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1846.208       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  3.95E-01  8.44E-01  1.37E+00  6.79E-01  9.27E-01  2.24E+00  1.00E-02  6.20E-01  8.87E-01  2.10E+00
 


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
+        1.17E+03
 
 TH 2
+       -4.00E+01  4.77E+02
 
 TH 3
+        7.34E+00  1.32E+02  5.63E+02
 
 TH 4
+       -2.66E+01  4.72E+02 -2.28E+02  9.96E+02
 
 TH 5
+        1.29E+01 -3.77E+02 -8.68E+02  1.02E+02  1.55E+03
 
 TH 6
+        1.84E+00 -7.39E+00  3.37E+00 -1.01E+01 -8.68E-01  2.16E+02
 
 TH 7
+        1.56E+00  3.37E+01  4.17E+00 -3.87E+00 -6.41E+00  4.53E-01  9.63E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.67E+00 -8.64E+00  7.02E-01 -4.16E+01  4.43E+01  1.89E+00  1.07E+01  0.00E+00  2.29E+02
 
 TH10
+        4.85E-01  1.18E+01 -2.67E+01 -1.07E+01 -4.04E+01  7.84E-01  1.45E+00  0.00E+00 -7.79E+00  7.50E+01
 
 TH11
+       -1.44E+01 -9.77E+00 -1.46E+01 -2.23E+01  5.14E+00  2.04E+00  1.52E+00  0.00E+00  1.43E+01  1.90E+01  1.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.393
Stop Time:
Wed Sep 29 21:54:28 CDT 2021
