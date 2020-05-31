# This dataset is created from the information in the file BHG_Seats_by_Section_05Sep19.xlsx

library(tidyverse)

tribble(
  ~location, ~proportion_of_stadium,
  1, 0.1866,
  2, 0.4825,
  3, 0.3189
) %>%
  write_tsv(., path = "stadium_seating.txt")
