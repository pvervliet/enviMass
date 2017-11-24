	
observeEvent(input$IS_drt2, {
    feedback(
      inputId = "IS_drt2",
      condition = input$IS_drt2 > 10,
	  text = "For most applications, this tolerance would be too large",
	  color = "#f23202"
    )
})

observeEvent(input$tar_drt2, {
    feedback(
      inputId = "tar_drt2",
      condition = input$tar_drt2 > 10,
	  text = "For most applications, this tolerance would be too large",
	  color = "#f23202"
    )
})

observeEvent(input$PWpath, {
	PWpath_ok <- TRUE
	if(
		!file.exists(input$PWpath)
	){
		PWpath_ok <- FALSE	
	}
    feedback(
      inputId = "PWpath",
      condition = !PWpath_ok,
	  text = "No msconvert.exe found for the specified path - please correct if you want to use msconvert for direct Thermo .raw file upload.",
	  color = "#f23202"
    )
})

  
	
