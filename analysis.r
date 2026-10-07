read_data <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d <- d[, c("Azimuth", "Lat", "Lon", "GpsLat", "GpsLon")]
  d <- d[!is.na(d$Azimuth), ]
  d$punto <- seq_len(nrow(d))
  d
}

d <- read_data("resection_log_golden.csv")


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
    r <- rbind(r, row_resecton(paste0("+0.05 lat/lon in P", i), dpos = d))
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
