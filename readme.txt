
In this short file, I help you to navigate within the code I provided. There are two "main" files, one for the central problem (questions 1 to 5 and 8) and one for the promoter characterization sub-problem (questions 6 and 7).
I use Boolean flags to activate or inactivate sections that are not relevant for some questions.

To answer question 2, 
1) you will need to edit the you_ode file that models the original (or reference) system. 
2) and you will edit the function main_dyn_range_optim that contains the code needed to plot the behaviors associated with the reference and extended models and ultimately find the circuit that has the largest dynamic range (hence the name).
You will set analyze_original_system to 1 and analyze_modified_systems to zero to work only with the reference model. You will edit only the first part of the analyze_original_system section

To answer question 3, edit the second part of the analyze_original_system section

To answer question 5, 
1) you will need to edit the you_odeR, you_odeI and you_odeRI files that corresponds to the  3 circuits we consider 
2) you will set analyze_modified_systems to 1 (and analyze_original_system to 0) and set use_default_params to 1 (and use_optimized_params to 0).

To answer question 6, go to main_prom_char, a function dedicated to the characterization of the promoter. Set select_char_circuit to 1 and infer_param_values to 0

To answer question 7, you ask me about the data. Set infer_param_values to 1 (and select_char_circuit to 0) in main_prom_char. 

To answer question 8, you go back to the main_dyn_range_optim function and work with analyze_modified_systems and use_optimized_params set to 1 (and use_default_params to 0)