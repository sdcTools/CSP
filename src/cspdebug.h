void   write_sol(char *);
void   write_heu(const char *);
int    control_ind(void);
void   print_rows(void);
int    print_row(CONSTRAINT *);
int    print_col(VARIABLE *);
void   print_sol();
int    control_constraint(CONSTRAINT *);
void   control_pool(void);
double violated_by_heur(CONSTRAINT *);
void   writing_lp(void);
void   print_card_cover(void);
void   testing_pool(void);

