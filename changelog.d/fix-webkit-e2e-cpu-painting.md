- Web tests: WebKit paints on the CPU (`WEBKIT_SKIA_ENABLE_CPU_RENDERING=1`). With GPU painting, the web
  process of a closed page died of SIGSEGV on a SkiaGPUWorker thread after its test had passed (5 times in
  2 local runs of the WebKit tests with the engine in the build, 0 in 3 runs on the CPU), and on WSL each
  crash left a dump of several GB on C:. The application is unchanged.
