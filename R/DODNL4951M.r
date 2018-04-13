DODNL4951M <- function(session){
updateSelectInput(session,"Ion_mode_Cal",selected = "none")
updateSelectInput(session,inputId="Cal_file_set",selected = "none")
updateNumericInput(session,"sel_meas_ID", value = 0)
updateNumericInput(session,"sel_meas", value = 0)
updateSelectInput(session,"Pos_compound_select",selected = "Choose")
updateSelectInput(session,"Neg_compound_select",selected = "Choose")
updateNumericInput(session,"sel_meas_comp", value = 0)
updateNumericInput(session,"sel_meas_comp_peak", value = 0)
updateNumericInput(session,"sel_meas_comp_comp", value = 0)
updateNumericInput(session,"atom_bound_peak", value = 0)
updateNumericInput(session,"sel_scans_ID", value = 0)
}
