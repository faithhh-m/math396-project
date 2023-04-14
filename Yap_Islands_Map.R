# Load necessary libraries
library(ggplot2)
library(mapdata)

# Get map data for Yap Islands
map <- map_data("worldHires")
yap_map <- subset(map, region == "Micronesia" & subregion == "Yap")

# Create a data frame with label locations and names
# ref: https://www.nejm.org/doi/full/10.1056/nejmoa0805715
label_data <- data.frame(
  long = c(138.17, 138.20, 138.17, 138.14, 138.115, 138.08, 138.10, 138.071352778, 138.07),
  lati = c(9.51, 9.53, 9.57, 9.56, 9.51, 9.475, 9.465, 9.44, 9.415),
  name = c("Tomil (21.5)", "Gagil (16.3)", "Maap (11.8)", "Fanif (8.4)", "Weloy (16.7)", "Dalipebinaw (6.2)", "Rull (11.2)", "Kanifay (3.6)", "Gilman (3.6)")
)

# Plot map using ggplot2
ggplot(yap_map, aes(x = long, y = lat)) +
  geom_polygon(fill = "lightblue", colour = "black") +
  geom_text(data = label_data, aes(x = long, y = lati, label = name), size = 4, fontface = "bold") +
  coord_fixed() +
  labs(title = "Attack Rates for Confirmed and Probable Zika Virus Disease per 1000",
       subtitle = "According to Municipality during the Period from April through July 2007") 
