---
source: /Users/byron/projects/rockbeatspaper/src/components/Card.tsx
relative: rockbeatspaper/src/components/Card.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { ReactNode } from "react";
import { cn } from "~/utils/cn";

interface CardProps {
  children: ReactNode;
  className?: string;
  header?: ReactNode;
  footer?: ReactNode;
}

export function Card({ children, className, header, footer }: CardProps) {
  return (
    <div
      className={cn(
        "bg-white dark:bg-gray-800 rounded-xl shadow-lg overflow-hidden",
        className
      )}
    >
      {header && (
        <div className="px-6 py-4 border-b border-gray-200 dark:border-gray-700">
          {header}
        </div>
      )}
      <div className="p-6">{children}</div>
      {footer && (
        <div className="px-6 py-4 border-t border-gray-200 dark:border-gray-700 bg-gray-50 dark:bg-gray-900">
          {footer}
        </div>
      )}
    </div>
  );
}

```
