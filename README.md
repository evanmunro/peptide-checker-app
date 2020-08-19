This is the source code for a [web application](https://evanmunro.shinyapps.io/peptide-checker-app/) that wraps a [C++ program](https://github.com/evanmunro/peptide-checker) written by the author for troubleshooting solid state peptide synthesis, so that those without C++ or command line knowledge can use the code. It contains one extra feature from the C++ utility, which is the ability to specify position specific protecting groups that adjust the mass of specific amino acids. 

Unfortunately, the free version of the shinyapps.io hosting that I'm using for web app has not been very reliable. As a result, I've also included instructions on how to run the application locally. 

### Running the Peptide Mass Checker Locally 

1. Make sure you have a current R installation, [instructions here](https://rstudio-education.github.io/hopr/starting.html) 
2. Run R from the command line

3. Install the R package `shiny` from R: 
```install.packages("shiny")```
4. Run the application

```shiny::runGitHub('peptide-checker-app', 'evanmunro')``` 
