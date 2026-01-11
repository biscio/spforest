#' Cross product for 3D vector
#'
#' @param a,b Vector c(x,y,z)
#'
#' @returns A vector with three coordinates. 
#' The function is different from the base R function crossprod 
#' which returns the scalar product between two vectors.
#' @export
#'
#' @examples
#' u = c(1, 2, 3)
#' v = c(2, 3, 4)
#' w = cross(u, v)
#' t(v)%*%w 
#' t(u)%*%w
cross <- function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' Random barycentric coordinates in a triangle
#'
#' @returns A vector with three coordinates.
#' @export
#'
#' @examples
#' random_barycentric()
random_barycentric <- function() {
  u <- runif(1)
  v <- runif(1)
  if (u + v > 1) {
    u <- 1 - u
    v <- 1 - v
  }
  c(u, v, 1 - u - v)
}


#' Build a 3D mesh from points and elevation data
#'
#' @param X A point pattern (of class "ppp")
#' @param elev An image (of class "im") representing the elevation
#' @param correction A numeric value to scale the elevation data for visual representation.
#'
#' @returns A list with two elements: \code{mesh}, a 3D mesh object, and \code{pp},
#' a matrix of 3D coordinates of the points projected onto the mesh.
#' @export
#'
#' @examples
#' res <- pptomesh(
#'   X = spatstat.data::bei,
#'   elev = spatstat.data::bei.extra$elev
#' )
pptomesh <- function(X, elev, correction = 6) {
  # TODO: Check that window of elev same as X. If note, what to do ?
  z <- correction * spatstat.geom::as.matrix.im(elev)
  nx <- ncol(z)
  ny <- nrow(z)

  vb <- t(cbind(
    expand.grid(
      x = elev$xcol,
      y = elev$yrow
    ),
    as.vector(t(z)),
    1
  ))

  # triangles
  faces <- matrix(0, nrow = 3, ncol = 2 * (nx - 1) * (ny - 1))
  k <- 1
  for (i in 1:(nx - 1)) {
    for (j in 1:(ny - 1)) {
      v0 <- (j - 1) * nx + i
      faces[, k] <- v0 + c(0, 1, nx)
      faces[, k + 1] <- v0 + c(1, nx + 1, nx)
      k <- k + 2
    }
  }

  mesh <- rgl::mesh3d(vertices = vb, triangles = faces)

  surface <- list(x = elev$xcol, y = elev$yrow, z = t(z))
  P <- matrix(c(X$x, X$y), ncol = 2)
  zP <- fields::interp.surface(surface, P)

  pointsech <- matrix(c(X$x, X$y, zP), ncol = 3)

  return(list(mesh = mesh, pp = pointsech))
}

#' Tirage des points de l'échantillon sur un mesh
#'
#' @param mesh A 3D mesh object (of class "mesh3d")
#' @param n Number of points to sample
#' @param weights Logical, whether to use weights for sampling
#' @param dimweight Dimension on which the weight are applied to (1 for x, 2 for y, 3 for z)
#' @param thres threshold value for weighting
#' @param p weight parameter
#'
#' @returns A list with two elements: \code{mesh}, the input mesh, and \code{pp}.
#' @export
#'
#' @examples
#' res <- pptomesh(
#'   X = spatstat.data::bei,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' dummypponmesh(
#'   mesh = res$mesh,
#'   n = 100,
#'   weights = TRUE,
#'   dimweight = 3,
#'   thres = 15
#' )
dummypponmesh <- function(mesh,
                          n,
                          weights = TRUE,
                          dimweight = 1,
                          thres = 15,
                          p = 0.7) {
  features <- features_mesh(mesh)
  tricenter <- features$tricenter
  # triarea <- features$triarea
  vertices <- features$vertices
  faces <- features$faces

  # Pondération : plus la coordonnée dimweight est grande, plus la probabilité est élevée
  if (weights == TRUE) {
    x_weights <- tricenter[, dimweight] # coordonnées x
    f <- function(u) {
      if (u >= thres) {
        return(p)
      }
      if (u < thres) {
        return(1 - p)
      }
    }
    x_weights <- sapply(x_weights, f)
    # x_weights <- x_weights - min(x_weights) # mettre à zéro le minimum
    # x_weights <- x_weights^2 # accentuer les grandes x, éviter 0
    # x_weights <- x_weights / max(x_weights)
    # Échantillonnage des triangles pondéré par x
    sampled_faces <- sample(1:nrow(faces),
      size = n,
      replace = TRUE,
      prob = x_weights
    )
  } else {
    sampled_faces <- sample(1:nrow(faces),
      size = n,
      replace = TRUE
    )
  }

  sampled_points <- t(sapply(sampled_faces, function(i) {
    A <- vertices[faces[i, 1], ]
    B <- vertices[faces[i, 2], ]
    C <- vertices[faces[i, 3], ]
    lambda <- random_barycentric()
    lambda[1] * A + lambda[2] * B + lambda[3] * C
  }))

  return(list(mesh = mesh, pp = sampled_points))
}

#' Tirage des points de l'échantillon sur un mesh (intern)
#'
#' @param tricenter Center of each triangle of the mesh
#' @param vertices Coordinates of each triange vertex in the mesh
#' @param faces Indices of each triangle vertex in the mesh
#' @param n Number of points to sample
#' @param weights Logical, whether to use weights for sampling
#' @param dimweight Dimension on which the weight are applied to (1 for x, 2 for y, 3 for z)
#'
#' @returns A matrix of sampled 3D points of size n x 3.
#' @export
#'
#' @examples
#' library(Rvcg)
#' data("humface")
#'
#' features <- features_mesh(humface)
#' tricenter <- features$tricenter
#' vertices <- features$vertices
#' faces <- features$faces
#' sample_points_manifold(tricenter,
#'   vertices,
#'   faces,
#'   n = 100,
#'   weights = TRUE,
#'   dimweight = 3
#' )
sample_points_manifold <- function(tricenter,
                                   vertices,
                                   faces,
                                   n,
                                   weights = TRUE,
                                   dimweight = 3) {
  # Pondération : plus la coordonnée dimweight est grande, plus la probabilité est élevée
  if (weights == TRUE) {
    x_weights <- tricenter[, dimweight] # coordonnées x
    x_weights <- x_weights - min(x_weights) # mettre à zéro le minimum
    x_weights <- x_weights^2 # accentuer les grandes x, éviter 0
    x_weights <- x_weights / max(x_weights)
    # Échantillonnage des triangles pondéré par x
    sampled_faces <- sample(1:nrow(faces), size = n, replace = TRUE, prob = x_weights)
  } else {
    sampled_faces <- sample(1:nrow(faces), size = n, replace = TRUE)
  }

  # Coordonnées barycentriques aléatoires dans chaque triangle
  random_barycentric <- function() {
    u <- runif(1)
    v <- runif(1)
    if (u + v > 1) {
      u <- 1 - u
      v <- 1 - v
    }
    c(u, v, 1 - u - v)
  }

  sampled_points <- t(sapply(sampled_faces, function(i) {
    A <- vertices[faces[i, 1], ]
    B <- vertices[faces[i, 2], ]
    C <- vertices[faces[i, 3], ]
    lambda <- random_barycentric()
    lambda[1] * A + lambda[2] * B + lambda[3] * C
  }))

  return(sampled_points)
}

#' Compute centres, areas, vertices, faces of triangular 3D mesh
#'
#' @param mesh A 3D mesh object (of class "mesh3d")
#'
#' @returns A list with four elements: \code{tricenter}, the centres of each triangle;
#' \code{triarea}, the area of each triangle; \code{vertices}, the coordinates of each vertex;
#' \code{faces}, the indices of each triangle vertex.
#' @export
#'
#' @examples
#' #' library(Rvcg)
#' data("humface")
#'
#' features <- features_mesh(humface)
features_mesh <- function(mesh) {
  # 1. Extraire les sommets du maillage
  vertices <- t(mesh$vb[1:3, ]) # N x 3
  faces <- t(mesh$it) # M x 3

  # 2. Calcul du centre de chaque triangle
  triangle_centers <- t(apply(faces, 1, function(f) {
    colMeans(vertices[f, ])
  }))

  # Calcul des aires de chaque triangle 
  # triangle_area <- function(A, B, C) {
  #   0.5 * norm(crossprod(matrix(B - A, ncol = 3), matrix(C - A, ncol = 3)), type = "2")
  # }
  # triangle_areas <- apply(faces, 1, function(f) {
  #   triangle_area(vertices[f[1], ], vertices[f[2], ], vertices[f[3], ])
  # })
  triangle_areas <- Rvcg::vcgArea(mesh, perface = TRUE)[[2]]
  return(list(
    tricenter = triangle_centers,
    triarea = triangle_areas,
    vertices = vertices,
    faces = faces
  ))
}

#' Tree intensity on 3D manifold
#'
#' @param intensity Number of dummy points to sample
#' @param triangle_centers Centres of each triangle of the mesh
#' @param triangle_areas Areas of each triangle of the mesh
#' @param vertices Coordinates of each triange vertex in the mesh
#' @param faces Indices of each triangle vertex in the mesh
#' @param pointsech Coordinates of each point of the observed point pattern on the mesh
#'
#' @returns A numeric vector of length equal to the number of triangles,
#' representing the estimated intensity for each triangle.
#' @export
#'
#' @examples
#' res <- pptomesh(
#'   X = spatstat.data::bei,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' features <- features_mesh(res$mesh)
#' manifold_tree(
#'   intensity = 200,
#'   tricenter = features$tricenter,
#'   triarea = features$triarea,
#'   vertices = features$vertices,
#'   faces = features$faces,
#'   pointsech = res$pp
#' )
manifold_tree <- function(intensity,
                          tricenter,
                          triarea,
                          vertices,
                          faces,
                          pointsech) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("The package RANN must be installed.")
  }

  # Echantillonage des dummy points
  points3D <- sample_points_manifold(
    tricenter,
    vertices,
    faces,
    n = intensity
  )
  # Trouver pour chaque triangle le point3D le plus proche
  nearest <- RANN::nn2(points3D, tricenter, k = 1) # k = 1 => plus proche

  # Indice du point3D le plus proche pour chaque triangle
  assigned_points <- nearest$nn.idx[, 1]

  # Regrouper les triangles par point
  # Liste : pour chaque point, les indices des triangles associés
  ntriangles <- nrow(tricenter)
  triangles_by_point <- split(1:ntriangles, assigned_points)

  # pour chaque dummy point, aire totale de sa partition
  partition_areas <- sapply(
    triangles_by_point,
    function(idx) sum(triarea[idx])
  )

  # Pour chaque dummy point, combien de points de ech sont les plus proches
  point_to_point <- RANN::nn2(points3D, pointsech, k = 1)$nn.idx[, 1]
  point_counts <- table(factor(point_to_point,
    levels = 1:nrow(points3D)
  ))

  # Calculer la densité pour chaque partition
  densities <- as.numeric(point_counts) / as.numeric(partition_areas)

  # Attribuer la densité à chaque triangle
  triangle_density <- numeric(ntriangles)
  for (i in seq_along(triangles_by_point)) {
    triangle_density[triangles_by_point[[i]]] <- densities[i]
  }

  return(triangle_density)
}

#' random forest intensity on mesh (intern)
#'
#' @param Ntrees Number of trees
#' @param intensity Number of dummy points to sample
#' @param mesh A 3D mesh object (of class "mesh3d")
#' @param pointsech Coordinates of each point of the observed point pattern on the mesh
#' @param verbose Logical. If TRUE show progress bar
#'
#' @returns A list of tree intensity estimates on manifold
#' @export
#'
#' @examples
#' res <- pptomesh(
#'   X = spatstat.data::bei,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' A <- manifold_forest(
#'   Ntrees = 10,
#'   intensity = 200,
#'   mesh = res$mesh,
#'   pointsech = res$pp,
#'   verbose = FALSE
#' )
manifold_forest <- function(Ntrees,
                            intensity,
                            mesh,
                            pointsech,
                            verbose = FALSE) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("The package rgl must be installed")
  }

  features <- features_mesh(mesh)
  tricenter <- features$tricenter
  triarea <- features$triarea
  vertices <- features$vertices
  faces <- features$faces

  # triangle_areas

  if (verbose) {
    progressr::handlers(global = TRUE)
    forestpgr <- function(x) {
      p <- progressr::progressor(along = x)
      lapply(1:Ntrees, FUN = function(j) {
        output <- manifold_tree(
          intensity = intensity,
          tricenter = tricenter,
          triarea = triarea,
          vertices = vertices,
          faces = faces,
          pointsech = pointsech
        )
        p(sprintf("x=%g", x))
        return(output)
      })
    }
    listmeshtree <- forestpgr(x = 1:Ntrees)
    progressr::handlers(global = FALSE)
  } else {
    listmeshtree <- lapply(1:Ntrees, FUN = function(j) {
      output <- manifold_tree(
        intensity = intensity,
        tricenter = tricenter,
        triarea = triarea,
        vertices = vertices,
        faces = faces,
        pointsech = pointsech
      )
      return(output)
    })
  }

  return(Reduce("+", listmeshtree) / Ntrees)

  # triangle_density <- manifold_tree(
  #   intensity = intensity,
  #   tricenter = tricenter,
  #   triarea = triarea,
  #   vertices = vertices,
  #   faces = faces,
  #   pointsech = pointsech
  # )
  # if (Ntrees > 1) {
  #   for (i in 2:Ntrees) {
  #     triangle_density <- triangle_density + manifold_tree(
  #       intensity = intensity,
  #       tricenter = tricenter,
  #       triarea = triarea,
  #       vertices = vertices,
  #       faces = faces,
  #       pointsech = pointsech
  #     )
  #   }
  # }
  #
  # return(triangle_density / Ntrees)
}


#' Plot manifold intensity (intern)
#'
#' @param forestmesh A spforestmesh object
#' @param points If TRUE, add the observed points on the plot
#' @param size description
#' @param colorbar If TRUE, add a colorbar to the plot
#' @param theta Angle for visualisation
#' @param phi Angle for visualisation
#' @param zoom Zoom for visualisation
#' @param pos Center position of the colorbar
#' @param nticks Number of ticks on the colorbar
#' @param lasttick If TRUE, add the last tick on the colorbar
#'
#' @returns A 3D plot, in a RGL window, of the manifold with intensity coloring
#' @export
#'
#' @examples
#' res <- pptomesh(
#'   X = spatstat.data::bei, ,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' forest <- spforest(X = res)
#' plot_manifold_intensity(log(forest + exp(-8)),
#'   points = TRUE,
#'   colorbar = TRUE,
#'   zoom = 0.65,
#'   nticks = 10,
#'   lasttick = TRUE,
#' )
plot_manifold_intensity <- function(forestmesh,
                                    points = FALSE,
                                    size = 2,
                                    colorbar = FALSE,
                                    theta = 10,
                                    phi = -50,
                                    zoom = 0.8,
                                    pos = NULL,
                                    nticks = 5,
                                    lasttick = TRUE) {
  triangle_density <- forestmesh$tridensity
  mesh <- forestmesh$mesh

  # Normaliser les densités entre 0 et 1
  normalized_density <- (triangle_density - min(triangle_density)) /
    (max(triangle_density) - min(triangle_density))

  # Obtenir la fonction de couleurs spatstat
  colfun <- spatstat.options("image.colfun") # c'est une fonction
  # Générer une palette avec, disons, 100 couleurs
  spatstat_palette <- colfun(100) # renvoie un vecteur de couleurs
  # Créer la fonction de dégradé pour interpolation continue
  col_fun <- colorRamp(spatstat_palette)
  # Appliquer aux densités normalisées
  rgb_colors <- col_fun(normalized_density) / 255

  # Créer un mesh avec couleurs par face
  colored_mesh <- rgl::tmesh3d(
    vertices = mesh$vb,
    indices = mesh$it,
    homogeneous = TRUE,
    material = list(color = rgb(rgb_colors)),
    meshColor = "faces"
  )

  # Affichage
  rgl::open3d(windowRect = c(100, 100, 900, 700)) # position + taille
  rgl::shade3d(colored_mesh)

  if (colorbar) {
    if (is.null(theta) || is.null(phi)) {
      stop("theta and phi must be provided")
    }
    rgl::view3d(theta = theta, phi = phi, zoom = zoom)
    # Centre de la colorbar (ajustable)
    vb <- t(colored_mesh$vb[1:3, ])
    if (is.null(pos)) {
      center <- c(max(vb[, 1]) * 1.1, min(vb[, 2]) +
        0.7 * diff(range(vb[, 2])), min(vb[, 3]))
    }
    add_colorbar3d(center,
      length = 3 / 5 * diff(range(vb[, 2])),
      width = 0.05 * diff(range(vb[, 1])),
      theta = theta,
      phi = phi,
      zlim = range(triangle_density),
      title = "",
      nticks = nticks,
      lasttick = lasttick,
      colfun = colfun
    )
  }

  if (points) {
    rgl::points3d(forestmesh$pp, col = "black", size = size, add = T)
  }
}


#' Create colorbar for 3D plot
#'
#' @param center Center position of the colorbar
#' @param length Size of color rectangles
#' @param width Width of the colorbar
#' @param theta Angle for visualisation
#' @param phi Angle for visualisation
#' @param zlim Limits of the colorbar
#' @param ncolors Number of colors in the colorbar
#' @param colfun Color function to generate colors
#' @param nticks Number of ticks on the colorbar
#' @param lasttick If TRUE, add the last tick on the colorbar
#' @param title Title of the colorbar
#' @param cex Character expansion for text
#'
#' @returns A colorbar added to the current RGL 3D plot
#' @export
add_colorbar3d <- function(center,
                           length = 1, width = 0.05,
                           theta = 0, phi = 0,
                           zlim = c(0, 1), ncolors = 100,
                           colfun = colorRampPalette(c("blue", "red")),
                           nticks, lasttick = TRUE,
                           title = NULL, cex = 1) {


  # Vecteurs caméra (dans l’espace)
  theta <- theta * pi / 180
  phi <- phi * pi / 180

  look <- c(cos(phi) * sin(theta), sin(phi), cos(phi) * cos(theta))
  up0 <- c(0, 1, 0)
  right <- cross(up0, look)
  right <- right / sqrt(sum(right^2))
  up <- cross(look, right)
  up <- up / sqrt(sum(up^2))

  # Couleurs
  colorscale <- colfun(ncolors)
  z_rel <- seq(-0.5, 0.5, length.out = ncolors + 1)

  # Dessiner les rectangles colorés
  for (i in 1:ncolors) {
    p1 <- center + z_rel[i] * length * up - (width / 2) * right
    p2 <- center + z_rel[i] * length * up + (width / 2) * right
    p3 <- center + z_rel[i + 1] * length * up + (width / 2) * right
    p4 <- center + z_rel[i + 1] * length * up - (width / 2) * right

    rgl::quads3d(rbind(p1, p2, p3, p4), col = colorscale[i], lit = FALSE)
  }

  # Ticks
  ticks <- pretty(zlim, n = nticks)
  if (lasttick) {
    nn <- length(ticks) + 1
  } else {
    nn <- length(ticks)
  }
  for (tick in ticks[-nn]) {
    zpos <- (tick - mean(zlim)) / diff(zlim) * length
    tpos <- center + zpos * up + 0.6 * width * right
    rgl::text3d(tpos[1], tpos[2], tpos[3],
      texts = format(tick, digits = 3),
      adj = c(0, 0.5), cex = cex
    )
  }

  # Titre
  if (!is.null(title)) {
    title_pos <- center + (length / 2 + 0.1 * length) * up
    rgl::text3d(title_pos[1], title_pos[2], title_pos[3],
      texts = title, adj = c(0.5, 0), cex = cex + 0.2
    )
  }
}
