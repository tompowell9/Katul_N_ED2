	  j,     k820309              13.0        éyT                                                                                                           
       farq_leuning.f90 FARQ_LEUNING                                                  
                   
                  ÿÿÿÿÿÿï                                                         HUGE                                                  
                 
                 »½×Ùß|Û=        1.D-10                                                                                    @              320#         @                                                     #LPHYSIOL_FULL%STOMA_DATA    #LPHYSIOL_FULL%DBLE    #CAN_PRSS    #CAN_RHOS    #CAN_SHV    #CAN_CO2    #IPFT    #LEAF_PAR    #LEAF_TEMP    #LINT_SHV    #GREEN_LEAF_FACTOR     #LEAF_AGING_FACTOR !   #LLSPAN "   #VM_BAR #   #LEAF_GBW $   #A_OPEN %   #A_CLOSED &   #GSW_OPEN '   #GSW_CLOSED (   #LSFC_SHV_OPEN )   #LSFC_SHV_CLOSED *   #LSFC_CO2_OPEN +   #LSFC_CO2_CLOSED ,   #LINT_CO2_OPEN -   #LINT_CO2_CLOSED .   #LEAF_RESP /   #VMOUT 0   #COMPPOUT 1   #LIMIT_FLAG 2   #OLD_ST_DATA 3                                                                                                                                                                                                                                                                                                                                   !@                               '@                    #RECALC    #T_L    #E_A 	   #PAR 
   #RB_FACTOR    #PRSS    #PHENOLOGY_FACTOR    #GSW_OPEN    #ILIMIT    #T_L_RESIDUAL    #E_A_RESIDUAL    #PAR_RESIDUAL    #RB_RESIDUAL    #PRSS_RESIDUAL    #LEAF_RESIDUAL    #GSW_RESIDUAL                                                                                                                                              1                                                              	                                              	               	                                              
               	                                                             	                                                             	                                                             	                                                             	                                                            	                                                      $       
   	                                                   (          	                                                   ,          	                                                   0          	                                                   4          	                                                   8          	                                                   <          	                                                   DBLE           
  @                                   	                
  @                                   	                
  @                                   	                
  @                                   	                
  @                                                    
  @                                   	                
  @                                   	                
  @                                   	                
  @                                    	                
  @                              !     	                
                                 "     	                
  @                              #     	                
  @                              $     	                D                                %     	                 D                                &     	                 D                                '     	                 D                                (     	                 D                                )     	                 D                                *     	                 D                                +     	                 D                                ,     	                 D                                -     	                 D                                .     	                 D                                /     	                 D                                0     	                 D                                1     	                 D @                               2                      
                                 3     @               #LPHYSIOL_FULL%STOMA_DATA    #         @                                  4                   #COMP_PHOTO_TEMPFUN%EXP 5   #COMP_PHOTO_TEMPFUN%MIN 6   #COMP_PHOTO_TEMPFUN%MAX 7   #COMP_PHOTO_TEMPFUN%DBLE 8   #IPFT 9   #LEAF_AGING_FACTOR :   #GREEN_LEAF_FACTOR ;                                                                                                                                                                               5     EXP                                            6     MIN                                            7     MAX                                            8     DBLE           
                                  9                     
  @                              :     	                
  @                              ;     	      #         @                                  <                    #IPFT =   #LIMIT_FLAG >             
                                  =                     D                                 >            %         @                               ?                   
       #ARRHENIUS%EXP @   #TEMP A   #REFVAL B   #HOR C                                                                        @     EXP           
                                 A     
                
                                 B     
                
                                 C     
      %         @                               D                    
       #TEMP E   #REFVAL F   #Q10 G             
                                 E     
                
                                 F     
                
                                 G     
      #         @                                  H                   #SET_CO2_DEMAND_PARAMS%TRIM I   #WHICHLIM J                                              I     TRIM           
  @                              J                    1 #         @                                  K                  #SOLVE_AOFIXED_CASE%SOLUTION_VARS L   #SOLVE_AOFIXED_CASE%SQRT S   #SOLVE_AOFIXED_CASE%MIN T   #ANSWER U   #SUCCESS V                     !@                           L     '0                    #LSFC_SHV M   #LSFC_CO2 N   #LINT_CO2 O   #CO2_DEMAND P   #STOM_COND_H2O Q   #STOM_COND_CO2 R                                              M                
                                              N               
                                              O               
                                              P               
                                              Q                
                                              R     (          
                                              S     SQRT                                            T     MIN           D                                 U     0               #SOLVE_AOFIXED_CASE%SOLUTION_VARS L             D                                 V            %         @                               W                     
       #         @                                  X                  #SOLVE_ITERATIVE_CASE%SOLUTION_VARS Y   #SOLVE_ITERATIVE_CASE%ABS `   #SOLVE_ITERATIVE_CASE%SQRT a   #SOLVE_ITERATIVE_CASE%MIN b   #SOLVE_ITERATIVE_CASE%DBLE c   #ANSWER d   #CONVERGED e                     !@                           Y     '0                    #LSFC_SHV Z   #LSFC_CO2 [   #LINT_CO2 \   #CO2_DEMAND ]   #STOM_COND_H2O ^   #STOM_COND_CO2 _                                              Z                
                                              [               
                                              \               
                                              ]               
                                              ^                
                                              _     (          
                                              `     ABS                                            a     SQRT                                            b     MIN                                            c     DBLE           D                                 d     0               #SOLVE_ITERATIVE_CASE%SOLUTION_VARS Y             D                                 e            %         @                               f                    
       #LINT_CO2 g             
                                 g     
      #         @                                  h                   #FIND_LINT_CO2_BOUNDS%SQRT i   #FIND_LINT_CO2_BOUNDS%MIN j   #FIND_LINT_CO2_BOUNDS%MAX k   #CIMIN l   #CIMAX m   #BOUNDED n                                              i     SQRT                                            j     MIN                                            k     MAX           D                                l     
                 D                                m     
                 D                                 n            #         @                                  o                    #NEWTON p   #LINT_CO2 q   #FUN r   #DERIV s             
                                  p                     
  @                              q     
                D                                r     
                 D                                s     
       %         @                               t                    
       #LINT_CO2 u   #CO2_DEMAND v             
                                 u     
                
                                 v     
      %         @                               w                    
       #LINT_CO2 x   #CO2_DEMAND y             
                                 x     
                
                                 y     
      %         @                               z                    
       #LINT_CO2 {   #STOM_COND_H2O |   #CO2_DEMAND }   #CO2_DEMAND_PRIME ~             
                                 {     
                
                                 |     
                
                                 }     
                
                                 ~     
             &      fn#fn    Æ   p       DISCARD    6  =       HUGE    s  v       TOLERFL8    é  s       MAXFPOFL    \  d      LPHYSIOL_FULL 6   À  ?     LPHYSIOL_FULL%STOMA_DATA+C34CONSTANTS =   ÿ  ¥   a   LPHYSIOL_FULL%STOMA_DATA%RECALC+C34CONSTANTS :   ¤  H   a   LPHYSIOL_FULL%STOMA_DATA%T_L+C34CONSTANTS :   ì  H   a   LPHYSIOL_FULL%STOMA_DATA%E_A+C34CONSTANTS :   4  H   a   LPHYSIOL_FULL%STOMA_DATA%PAR+C34CONSTANTS @   |  H   a   LPHYSIOL_FULL%STOMA_DATA%RB_FACTOR+C34CONSTANTS ;   Ä  H   a   LPHYSIOL_FULL%STOMA_DATA%PRSS+C34CONSTANTS G   	  H   a   LPHYSIOL_FULL%STOMA_DATA%PHENOLOGY_FACTOR+C34CONSTANTS ?   T	  H   a   LPHYSIOL_FULL%STOMA_DATA%GSW_OPEN+C34CONSTANTS =   	  H   a   LPHYSIOL_FULL%STOMA_DATA%ILIMIT+C34CONSTANTS C   ä	  H   a   LPHYSIOL_FULL%STOMA_DATA%T_L_RESIDUAL+C34CONSTANTS C   ,
  H   a   LPHYSIOL_FULL%STOMA_DATA%E_A_RESIDUAL+C34CONSTANTS C   t
  H   a   LPHYSIOL_FULL%STOMA_DATA%PAR_RESIDUAL+C34CONSTANTS B   ¼
  H   a   LPHYSIOL_FULL%STOMA_DATA%RB_RESIDUAL+C34CONSTANTS D     H   a   LPHYSIOL_FULL%STOMA_DATA%PRSS_RESIDUAL+C34CONSTANTS D   L  H   a   LPHYSIOL_FULL%STOMA_DATA%LEAF_RESIDUAL+C34CONSTANTS C     H   a   LPHYSIOL_FULL%STOMA_DATA%GSW_RESIDUAL+C34CONSTANTS #   Ü  =      LPHYSIOL_FULL%DBLE '     @   a   LPHYSIOL_FULL%CAN_PRSS '   Y  @   a   LPHYSIOL_FULL%CAN_RHOS &     @   a   LPHYSIOL_FULL%CAN_SHV &   Ù  @   a   LPHYSIOL_FULL%CAN_CO2 #     @   a   LPHYSIOL_FULL%IPFT '   Y  @   a   LPHYSIOL_FULL%LEAF_PAR (     @   a   LPHYSIOL_FULL%LEAF_TEMP '   Ù  @   a   LPHYSIOL_FULL%LINT_SHV 0     @   a   LPHYSIOL_FULL%GREEN_LEAF_FACTOR 0   Y  @   a   LPHYSIOL_FULL%LEAF_AGING_FACTOR %     @   a   LPHYSIOL_FULL%LLSPAN %   Ù  @   a   LPHYSIOL_FULL%VM_BAR '     @   a   LPHYSIOL_FULL%LEAF_GBW %   Y  @   a   LPHYSIOL_FULL%A_OPEN '     @   a   LPHYSIOL_FULL%A_CLOSED '   Ù  @   a   LPHYSIOL_FULL%GSW_OPEN )     @   a   LPHYSIOL_FULL%GSW_CLOSED ,   Y  @   a   LPHYSIOL_FULL%LSFC_SHV_OPEN .     @   a   LPHYSIOL_FULL%LSFC_SHV_CLOSED ,   Ù  @   a   LPHYSIOL_FULL%LSFC_CO2_OPEN .     @   a   LPHYSIOL_FULL%LSFC_CO2_CLOSED ,   Y  @   a   LPHYSIOL_FULL%LINT_CO2_OPEN .     @   a   LPHYSIOL_FULL%LINT_CO2_CLOSED (   Ù  @   a   LPHYSIOL_FULL%LEAF_RESP $     @   a   LPHYSIOL_FULL%VMOUT '   Y  @   a   LPHYSIOL_FULL%COMPPOUT )     @   a   LPHYSIOL_FULL%LIMIT_FLAG *   Ù  f   a   LPHYSIOL_FULL%OLD_ST_DATA #   ?  r      COMP_PHOTO_TEMPFUN '   ±  <      COMP_PHOTO_TEMPFUN%EXP '   í  <      COMP_PHOTO_TEMPFUN%MIN '   )  <      COMP_PHOTO_TEMPFUN%MAX (   e  =      COMP_PHOTO_TEMPFUN%DBLE (   ¢  @   a   COMP_PHOTO_TEMPFUN%IPFT 5   â  @   a   COMP_PHOTO_TEMPFUN%LEAF_AGING_FACTOR 5   "  @   a   COMP_PHOTO_TEMPFUN%GREEN_LEAF_FACTOR ,   b  b       PHOTOSYNTHESIS_EXACT_SOLVER 1   Ä  @   a   PHOTOSYNTHESIS_EXACT_SOLVER%IPFT 7     @   a   PHOTOSYNTHESIS_EXACT_SOLVER%LIMIT_FLAG    D         ARRHENIUS    à  <      ARRHENIUS%EXP      @   a   ARRHENIUS%TEMP !   \  @   a   ARRHENIUS%REFVAL      @   a   ARRHENIUS%HOR    Ü  o       COLLATZ    K  @   a   COLLATZ%TEMP      @   a   COLLATZ%REFVAL    Ë  @   a   COLLATZ%Q10 &     v       SET_CO2_DEMAND_PARAMS +     =      SET_CO2_DEMAND_PARAMS%TRIM /   ¾  L   a   SET_CO2_DEMAND_PARAMS%WHICHLIM #   
  À       SOLVE_AOFIXED_CASE >   Ê  °      SOLVE_AOFIXED_CASE%SOLUTION_VARS+C34CONSTANTS G   z  H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%LSFC_SHV+C34CONSTANTS G   Â  H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%LSFC_CO2+C34CONSTANTS G   
  H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%LINT_CO2+C34CONSTANTS I   R  H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%CO2_DEMAND+C34CONSTANTS L     H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%STOM_COND_H2O+C34CONSTANTS L   â  H   a   SOLVE_AOFIXED_CASE%SOLUTION_VARS%STOM_COND_CO2+C34CONSTANTS (   *  =      SOLVE_AOFIXED_CASE%SQRT '   g  <      SOLVE_AOFIXED_CASE%MIN *   £  n   a   SOLVE_AOFIXED_CASE%ANSWER +     @   a   SOLVE_AOFIXED_CASE%SUCCESS "   Q  P       FIND_TWILIGHT_MIN %   ¡        SOLVE_ITERATIVE_CASE @   ¦   °      SOLVE_ITERATIVE_CASE%SOLUTION_VARS+C34CONSTANTS I   V!  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%LSFC_SHV+C34CONSTANTS I   !  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%LSFC_CO2+C34CONSTANTS I   æ!  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%LINT_CO2+C34CONSTANTS K   ."  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%CO2_DEMAND+C34CONSTANTS N   v"  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%STOM_COND_H2O+C34CONSTANTS N   ¾"  H   a   SOLVE_ITERATIVE_CASE%SOLUTION_VARS%STOM_COND_CO2+C34CONSTANTS )   #  <      SOLVE_ITERATIVE_CASE%ABS *   B#  =      SOLVE_ITERATIVE_CASE%SQRT )   #  <      SOLVE_ITERATIVE_CASE%MIN *   »#  =      SOLVE_ITERATIVE_CASE%DBLE ,   ø#  p   a   SOLVE_ITERATIVE_CASE%ANSWER /   h$  @   a   SOLVE_ITERATIVE_CASE%CONVERGED     ¨$  ^       CALC_CO2_DEMAND )   %  @   a   CALC_CO2_DEMAND%LINT_CO2 %   F%  Æ       FIND_LINT_CO2_BOUNDS *   &  =      FIND_LINT_CO2_BOUNDS%SQRT )   I&  <      FIND_LINT_CO2_BOUNDS%MIN )   &  <      FIND_LINT_CO2_BOUNDS%MAX +   Á&  @   a   FIND_LINT_CO2_BOUNDS%CIMIN +   '  @   a   FIND_LINT_CO2_BOUNDS%CIMAX -   A'  @   a   FIND_LINT_CO2_BOUNDS%BOUNDED !   '  v       ITER_SOLVER_STEP (   ÷'  @   a   ITER_SOLVER_STEP%NEWTON *   7(  @   a   ITER_SOLVER_STEP%LINT_CO2 %   w(  @   a   ITER_SOLVER_STEP%FUN '   ·(  @   a   ITER_SOLVER_STEP%DERIV #   ÷(  n       CALC_STOM_COND_H2O ,   e)  @   a   CALC_STOM_COND_H2O%LINT_CO2 .   ¥)  @   a   CALC_STOM_COND_H2O%CO2_DEMAND &   å)  n       CALC_CO2_DEMAND_PRIME /   S*  @   a   CALC_CO2_DEMAND_PRIME%LINT_CO2 1   *  @   a   CALC_CO2_DEMAND_PRIME%CO2_DEMAND )   Ó*         CALC_STOM_COND_H2O_PRIME 2   j+  @   a   CALC_STOM_COND_H2O_PRIME%LINT_CO2 7   ª+  @   a   CALC_STOM_COND_H2O_PRIME%STOM_COND_H2O 4   ê+  @   a   CALC_STOM_COND_H2O_PRIME%CO2_DEMAND :   *,  @   a   CALC_STOM_COND_H2O_PRIME%CO2_DEMAND_PRIME 