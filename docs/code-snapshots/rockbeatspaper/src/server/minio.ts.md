---
source: /Users/byron/projects/rockbeatspaper/src/server/minio.ts
relative: rockbeatspaper/src/server/minio.ts
generated_at: 2025-12-23 10:28
---

```typescript
import { Client } from "minio";
import { env } from "./env";
import { getBaseUrl } from "./utils/base-url";

export const minioBaseUrl = getBaseUrl({ port: 9000 });

export const minioClient = new Client({
  endPoint: minioBaseUrl.split("://")[1]!,
  useSSL: minioBaseUrl.startsWith("https://"),
  accessKey: "admin",
  secretKey: env.ADMIN_PASSWORD,
});

```
