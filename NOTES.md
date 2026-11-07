# Kurz-Notes zu A5

## Setup
- `app/searchSpeedTest.cpp` erstellt (kopiert von `main2.cpp`, umbenannt)
- `app/CMakeLists.txt` um Target `searchSpeedTest` erweitert:
  ```cmake
  add_executable(searchSpeedTest searchSpeedTest.cpp)
  add_dependencies(searchSpeedTest learnSPH)
  target_link_libraries(searchSpeedTest learnSPH)
  ```
- CMake neu konfigurieren nach CMakeLists-Änderung: `cmake ../framework` im Build-Verzeichnis

## Brute Force
- `learnSPH/bruteforcesearch.h` erweitert: `find_neighbors(const vector<Vector3d>&, double radius, vector<vector<unsigned int> >&)`
- Nutzt `d2 = dx*dx+dy*dy+dz*dz <= r*r` (vermeidet sqrt)
- Getestet mit 20x20x20 Gitter (8000 Partikel)
- Diskutiere ersten Run verwerfen, mehrere Runs mitteln

## CompactNSearch
- API: `NeighborhoodSearch nsearch(radius); nsearch.add_point_set(x.data(), n, true); nsearch.find_neighbors();`
- Mit 50x50x50 (125000) Partikel testen
- Nachbarn abfragen: `ps.n_neighbors(0, i)`

## Timing
- `std::chrono::high_resolution_clock`
- 1x ohne Messen, dann 3-5x mitteln
