# %% [markdown]

read_data <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[, c("Azimuth", "Lat", "Lon", "GpsLat", "GpsLon")]
  d <- d[!is.na(d$Azimuth), ]
  d$punto <- seq_len(nrow(d))
  d
}

d <- read_data("resection_log_golden1.csv")


resection <- function(az, x, y, add = rep(TRUE, length(az))) {
  az <- az[add]
  x <- x[add]
  y <- y[add]
  a <- az * rad
  nx <- cos(a)
  ny <- -sin(a)
  b <- nx * x + ny * y
  A <- matrix(c(sum(nx^2), sum(nx * ny), sum(nx * ny), sum(ny^2)), 2, 2)
  p <- solve(A, c(sum(nx * b), sum(ny * b)))
  rex <- nx * (p[1] - x) + ny * (p[2] - y)
  list(east = p[1], north = p[2], error = sqrt(sum(p^2)))
}


RADIUS_EARTH <- 6371008.0
rad <- pi / 180

to_degrees <- function(x, y, origin) {
  c(
    lat = origin[["lat"]] + y / (rad * RADIUS_EARTH),
    lon = origin[["lon"]] +
      x / (rad * RADIUS_EARTH * cos(origin[["lat"]] * rad))
  )
}

to_meters <- function(lat, lon, origin) {
  x <- (lon - origin[["lon"]]) * rad * RADIUS_EARTH * cos(origin[["lat"]] * rad)
  y <- (lat - origin[["lat"]]) * rad * RADIUS_EARTH
  cbind(x = x, y = y)
}

run <- function(table, origin, k = 3, daz_val = 1, dpos_val = 0.0001) {
  t <- table[seq_len(k), ]
  n <- nrow(t)

  row_resecton <- function(
    name,
    daz = rep(0, n),
    dpos = rep(0, n),
    add = rep(TRUE, n)
  ) {
    lat <- t$lat + dpos
    lon <- t$lon + dpos
    az <- t$az + daz
    xy <- to_meters(lat, lon, origin)
    s <- resection(az, xy[, "x"], xy[, "y"], add)
    sol <- to_degrees(s$east, s$north, origin)
    P <- ifelse(add, round(az, 2), NA)
    Q <- ifelse(add, sprintf("(%.6f, %.6f)", lat, lon), NA)
    out <- as.data.frame(
      c(
        list(scenario = name),
        setNames(as.list(P), paste0("P", 1:n)),
        setNames(as.list(Q), paste0("Q", 1:n))
      ),
      stringAsFactors = FALSE
    )

    out$Result1 <- sprintf("(%.7f, %.7f)", sol[["lat"]], sol[["lon"]])
    out$Result2 <- sprintf("(%.7f, %.7f)", origin[["lat"]], origin[["lon"]])
    out$error_origin_m <- s$error
    out
  }
  r <- row_resecton("BASE")
  for (i in 1:n) {
    d <- rep(0, n)
    d[i] <- daz_val
    r <- rbind(r, row_resecton(paste0("+1 Azimuth angle in P", i), daz = d))
  }
  r <- rbind(r, row_resecton("+1 Azimuth angle in all", daz = rep(daz_val, n)))
  for (i in 1:n) {
    d <- rep(0, n)
    d[i] <- dpos_val
    r <- rbind(r, row_resecton(paste0("+0.0001 lat/lon in P", i), dpos = d))
  }
  r <- rbind(r, row_resecton("+0.05 lat/lon in all", dpos = rep(dpos_val, n)))
  if (n > 2) {
    for (i in 1:n) {
      r <- rbind(
        r,
        row_resecton(paste0("Without P", i), add = seq_len(n) != i)
      )
    }
  }
  r
}

origin <- c(lat = 48.7022484, lon = 2.1530186)
table <- data.frame(point = d$punto, az = d$Azimuth, lat = d$Lat, lon = d$Lon)

res <- run(table, origin, k = 3)

print(d)
print(res)
View(res)

# ---- requiere en memoria: tabla (az, lat, lon), origen, res (tabla de 10 columnas), a_metros() ----
install.packages("ggplot2")
library(ggplot2)
Sys.setlocale("LC_ALL", "English_United States.UTF-8")
AZUL <- "#2a78d6"
NARANJA <- "#eb6834"
VERDE <- "#1baf7a"
GRIS <- "#8a8a8a"
dir.create("charts", showWarnings = FALSE)


# ---------- datos comunes ----------
xy <- to_meters(table$lat, table$lon, origin)
pines <- data.frame(
  punto = paste0("P", table$punto),
  x = xy[, "x"],
  y = xy[, "y"],
  az = table$az
)

parse_ll <- function(s) {
  n <- as.numeric(strsplit(gsub("[()]", "", s), ", ")[[1]])
  c(lat = n[1], lon = n[2])
}
sol_m <- t(sapply(res$Result1, function(s) {
  p <- parse_ll(s)
  to_meters(p[["lat"]], p[["lon"]], origin)
}))
dimnames(sol_m) <- list(NULL, c("x", "y"))
base <- sol_m[1, ]


# ---------- 1. MAPA EN METROS ----------
L <- 4000
rayos <- transform(pines, dx = sin(az * pi / 180), dy = cos(az * pi / 180))
rayos$x0 <- rayos$x - L * rayos$dx
rayos$y0 <- rayos$y - L * rayos$dy
rayos$x1 <- rayos$x + L * rayos$dx
rayos$y1 <- rayos$y + L * rayos$dy

todos_x <- c(pines$x, 0, base[["x"]])
todos_y <- c(pines$y, 0, base[["y"]])
mx <- 0.12 * diff(range(todos_x))
my <- 0.12 * diff(range(todos_y))

g1 <- ggplot() +
  geom_segment(
    data = rayos,
    aes(
      x = x0,
      y = y0,
      xend = x1,
      yend = y1,
      linetype = "Recta medida (azimut)",
      colour = "Recta medida (azimut)"
    ),
    linewidth = 0.6
  ) +
  geom_segment(
    data = pines,
    aes(
      x = 0,
      y = 0,
      xend = x,
      yend = y,
      linetype = "Dirección real (origen → pin)",
      colour = "Dirección real (origen → pin)"
    ),
    linewidth = 0.8
  ) +
  geom_point(data = pines, aes(x, y), size = 3) +
  geom_text(
    data = pines,
    aes(x, y, label = punto),
    nudge_y = 0.03 * diff(range(todos_y)),
    fontface = "bold",
    vjust = 0
  ) +
  geom_point(aes(x = 0, y = 0), shape = 8, size = 4) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = "origen real",
    hjust = 1.15,
    vjust = -0.9,
    size = 3.6
  ) +
  geom_point(aes(x = base[["x"]], y = base[["y"]]), colour = AZUL, size = 4) +
  annotate(
    "text",
    x = base[["x"]],
    y = base[["y"]],
    colour = AZUL,
    fontface = "bold",
    size = 3.8,
    label = sprintf("solución (BASE)\nerror %.0f m", res$error_origin_m[1]),
    hjust = -0.12,
    vjust = 0.5
  ) +
  scale_linetype_manual(
    values = c(
      "Recta medida (azimut)" = "solid",
      "Dirección real (origen → pin)" = "dashed"
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "Recta medida (azimut)" = GRIS,
      "Dirección real (origen → pin)" = VERDE
    ),
    name = NULL
  ) +
  coord_fixed(
    xlim = range(todos_x) + c(-mx, mx),
    ylim = range(todos_y) + c(-my, my)
  ) +
  labs(
    title = "Mapa en metros: rectas medidas, dirección real y solución",
    x = "este (m)",
    y = "norte (m)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave("charts/1_map_metros.png", g1, width = 8, height = 7, dpi = 150)
