---
source: /Users/byron/projects/rockbeatspaper/tsconfig.json
relative: rockbeatspaper/tsconfig.json
generated_at: 2025-12-23 10:28
---

```json
{
  "include": [
    "**/*.ts",
    "**/*.tsx",
    "**/*.mjs",
    "**/*.cjs",
    "**/*.js",
    "**/*.jsx"
  ],
  "compilerOptions": {
    "strict": true,
    "noUncheckedIndexedAccess": true,
    "esModuleInterop": true,
    "jsx": "react-jsx",
    "target": "ESNext",
    "moduleResolution": "Bundler",
    "skipLibCheck": true,
    "lib": ["DOM", "DOM.Iterable", "ES2022"],
    "allowJs": true,
    "forceConsistentCasingInFileNames": true,
    "baseUrl": ".",
    "paths": {
      "~/*": ["./src/*"]
    },
    "noEmit": true
  }
}

```
