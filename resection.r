data <- read.csv("resection_log_20260928_003444.csv", stringsAsFactors = FALSE)

str(data)

data <- data[, c("Azimuth", "Lat", "Lon", "GpsLat", "GpsLon")]

data


read_log <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[, c(
    "Timestamp",
    "SessionNumber",
    "Azimuth",
    "Lat",
    "Lon",
    "GpsLat",
    "GpsLon"
  )]
  d$file <- basename(path)
  d$session <- paste(d$file, d$SessionNumber, sep = "_S")
  d$date <- as.Date(substr(d$Timestamp, 1, 10))
  d
}

files <- list.files(pattern = "^resection_log_.*\\.csv$")
files2 <- list.files(pattern = "^resection_log_golden*\\.csv$")
data <- do.call(rbind, lapply(files2, read_log))


sessions <- split(data, data$session)
sessions

process_session <- function(d) {
  is_marker <- !is.na(d$Azimuth)
  is_pin <- !is.na(d$Lat) & is.na(d$GpsLat)

  pins_before <- cumsum(is_pin) - is_pin

  markers <- data.frame(
    point = pins_before[is_marker] + 1,
    az = d$Azimuth[is_marker]
  )

  pins <- data.frame(
    point = seq_len(sum(is_pin)),
    lat = d$Lat[is_pin],
    lon = d$Lon[is_pin]
  )

  table <- merge(markers, pins, by = "point")

  origin <- c(
    lat = mean(d$GpsLat[is_marker]),
    lon = mean(d$GpsLon[is_marker])
  )

  list(table = table, origin = origin)
}

r <- process_session(sessions[[1]])

table <- r$table
origin <- r$origin

metersPerDegLat <- 111320.0

RADIUS_EARTH <- 6371008.0
rad <- pi / 180

to_meters <- function(lat, lon, origin) {
  x <- (lon - origin[["lon"]]) * rad * RADIUS_EARTH * cos(origin[["lat"]] * rad)
  y <- (lat - origin[["lat"]]) * rad * RADIUS_EARTH
  cbind(x = x, y = y)
}


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

DAZ <- 1
DPOS <- 0.05

n <- nrow(table)
un <- function(k) c(seq_len(n) == k)

sensitivity_table <- function(table, origin, k, daz_val = 1, dpos_val = 0.00001) {
  t <- table[seq_len(k), ]
  n <- nrow(t)

  run_scenario <- function(
    name,
    daz = rep(0, n),
    dpos = rep(0, n),
    add = rep(TRUE, n)
  ) {
    xy <- to_meters(table$lat + dpos, table$lon + dpos, origin)
    s <- resection(table$az + daz, xy[, "x"], xy[, "y"], add)
    data.frame(
      scenario = name,
      east_m = s$east,
      north_m = s$north,
      error_origin_m = s$error
    )
  }

  sens <- run_scenario("BASE")

  for (k in 1:n) {
    sens <- rbind(
      sens,
      run_scenario(paste0("+1 Azimuth angle in P", k), daz = DAZ * un(k))
    )
  }
  sens <- rbind(
    sens,
    run_scenario("+1 Azimuth angle in all", daz = rep(DAZ, n))
  )
  for (k in 1:n) {
    sens <- rbind(
      sens,
      run_scenario(paste0("+0.05 lat/lon in P", k), dpos = DPOS * un(k))
    )
  }
  sens <- rbind(
    sens,
    run_scenario("+0.05 lat/lon in all", dpos = rep(DPOS, n))
  )

  # "Sin Pk" solo si quedan al menos 2 rayos
  if (n > 2) {
    for (i in 1:n) {
      sens <- rbind(
        sens,
        run_scenario(paste0("Without P", i), add = seq_len(n) != i)
      )
    }
  }

  sens$delta_m <- sqrt(
    (sens$east_m - sens$east_m[1])^2 + (sens$north_m - sens$north_m[1])^2
  )
  sens$points <- k

  sens
}

tables <- lapply(2:5, function(k) sensitivity_table(table, origin, k))
names(tables) <- paste0("k", 2:5)

tables$k2
View(tables$k5)
View(tables$k4)
View(tables$k3)
View(tables$k2)
View(sapply(tables, nrow))

tables2 <- lapply(2:5, function(k) sensitivity_table(table, origin, k))
names(tables2) <- paste0("k", 2:5)

View(tables2)


data
sens
View(sens)
