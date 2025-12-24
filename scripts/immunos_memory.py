#!/usr/bin/env python3
"""
IMMUNOS Memory System (T Cell Pattern)

Adaptive memory system for Claude Code persistence. Stores conversations,
decisions, and context with priority-based decay.

Usage:
    python3 immunos_memory.py create --content "..." --priority high
    python3 immunos_memory.py store --content "..." --priority high
    python3 immunos_memory.py list --priority high
    python3 immunos_memory.py search "keyword" --type conversation
    python3 immunos_memory.py export --type conversation --output memories.json
    python3 immunos_memory.py decay --dry-run
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime, timedelta
import hashlib


class Memory:
    """Single memory unit with adaptive decay"""

    PRIORITY_DECAY_RATES = {
        'high': 0.01,      # Slow decay - critical information
        'medium': 0.1,     # Moderate decay - useful context
        'low': 0.5,        # Fast decay - temporary info
    }

    PRIORITY_LEVELS = ['high', 'medium', 'low']

    def __init__(self,
                 content: Dict,
                 priority: str = 'medium',
                 memory_type: str = 'conversation',
                 references: Optional[List[str]] = None,
                 tags: Optional[List[str]] = None,
                 memory_id: Optional[str] = None,
                 timestamp: Optional[str] = None,
                 relevance_score: float = 1.0,
                 access_count: int = 0):

        self.content = content
        self.priority = priority if priority in self.PRIORITY_LEVELS else 'medium'
        self.memory_type = memory_type
        self.references = references or []
        self.tags = tags or []
        self.memory_id = memory_id or self._generate_id()
        self.timestamp = timestamp or datetime.now().isoformat()
        self.relevance_score = relevance_score
        self.access_count = access_count
        self.decay_rate = self.PRIORITY_DECAY_RATES[self.priority]
        self.last_accessed = self.timestamp

    def _generate_id(self) -> str:
        """Generate unique memory ID"""
        timestamp = datetime.now().isoformat()
        content_str = json.dumps(self.content, sort_keys=True)
        hash_input = f"{timestamp}{content_str}".encode()
        return f"mem_{hashlib.sha256(hash_input).hexdigest()[:16]}"

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization"""
        return {
            'memory_id': self.memory_id,
            'timestamp': self.timestamp,
            'type': self.memory_type,
            'priority': self.priority,
            'decay_rate': self.decay_rate,
            'relevance_score': self.relevance_score,
            'access_count': self.access_count,
            'last_accessed': self.last_accessed,
            'content': self.content,
            'references': self.references,
            'tags': self.tags,
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'Memory':
        """Create Memory from dictionary"""
        return cls(
            content=data['content'],
            priority=data['priority'],
            memory_type=data['type'],
            references=data.get('references', []),
            tags=data.get('tags', []),
            memory_id=data['memory_id'],
            timestamp=data['timestamp'],
            relevance_score=data.get('relevance_score', 1.0),
            access_count=data.get('access_count', 0),
        )

    def decay(self, hours_elapsed: float):
        """Apply decay to relevance score"""
        decay_amount = self.decay_rate * (hours_elapsed / 24.0)  # Decay per day
        self.relevance_score = max(0.0, self.relevance_score - decay_amount)

    def reinforce(self, boost: float = 0.1):
        """Reinforce memory (accessed or relevant)"""
        self.relevance_score = min(1.0, self.relevance_score + boost)
        self.access_count += 1
        self.last_accessed = datetime.now().isoformat()

    def is_expired(self, threshold: float = 0.1) -> bool:
        """Check if memory has decayed below threshold"""
        return self.relevance_score < threshold


class MemoryStore:
    """Storage and retrieval system for memories"""

    def __init__(self, base_path: str = '/Users/byron/projects/.immunos'):
        self.base_path = Path(base_path)
        self.memory_path = self.base_path / 'memory'
        self.conversations_path = self.memory_path / 'conversations'
        self.decisions_path = self.memory_path / 'decisions'
        self.preferences_path = self.memory_path / 'preferences'
        self.index_path = self.memory_path / 'index.json'

        # Ensure directories exist
        for path in [self.conversations_path, self.decisions_path, self.preferences_path]:
            path.mkdir(parents=True, exist_ok=True)

        # Load or create index
        self.index = self._load_index()

    def _load_index(self) -> Dict:
        """Load memory index"""
        if self.index_path.exists():
            with open(self.index_path, 'r') as f:
                return json.load(f)
        return {
            'memories': {},
            'by_priority': {'high': [], 'medium': [], 'low': []},
            'by_type': {},
            'by_tag': {},
            'stats': {
                'total_memories': 0,
                'active_memories': 0,
                'expired_memories': 0,
                'last_decay': None,
            }
        }

    def _save_index(self):
        """Save memory index"""
        with open(self.index_path, 'w') as f:
            json.dump(self.index, f, indent=2)

    def _get_storage_path(self, memory: Memory) -> Path:
        """Determine storage path based on memory type"""
        if memory.memory_type == 'decision':
            return self.decisions_path
        elif memory.memory_type == 'preference':
            return self.preferences_path
        else:
            return self.conversations_path

    def create_memory(self, memory: Memory) -> str:
        """Store a new memory"""
        # Save memory file
        storage_path = self._get_storage_path(memory)
        memory_file = storage_path / f"{memory.memory_id}.json"

        with open(memory_file, 'w') as f:
            json.dump(memory.to_dict(), f, indent=2)

        # Update index
        self.index['memories'][memory.memory_id] = {
            'type': memory.memory_type,
            'priority': memory.priority,
            'timestamp': memory.timestamp,
            'file': str(memory_file.relative_to(self.memory_path)),
            'relevance_score': memory.relevance_score,
            'tags': memory.tags,
        }

        # Update priority index
        if memory.memory_id not in self.index['by_priority'][memory.priority]:
            self.index['by_priority'][memory.priority].append(memory.memory_id)

        # Update type index
        if memory.memory_type not in self.index['by_type']:
            self.index['by_type'][memory.memory_type] = []
        if memory.memory_id not in self.index['by_type'][memory.memory_type]:
            self.index['by_type'][memory.memory_type].append(memory.memory_id)

        # Update tag index
        for tag in memory.tags:
            if tag not in self.index['by_tag']:
                self.index['by_tag'][tag] = []
            if memory.memory_id not in self.index['by_tag'][tag]:
                self.index['by_tag'][tag].append(memory.memory_id)

        # Update stats
        self.index['stats']['total_memories'] += 1
        self.index['stats']['active_memories'] += 1

        self._save_index()

        return memory.memory_id

    def get_memory(self, memory_id: str) -> Optional[Memory]:
        """Retrieve a specific memory"""
        if memory_id not in self.index['memories']:
            return None

        memory_info = self.index['memories'][memory_id]
        memory_file = self.memory_path / memory_info['file']

        if not memory_file.exists():
            return None

        with open(memory_file, 'r') as f:
            data = json.load(f)

        memory = Memory.from_dict(data)

        # Reinforce on access
        memory.reinforce(boost=0.05)
        self.update_memory(memory)

        return memory

    def update_memory(self, memory: Memory):
        """Update an existing memory"""
        if memory.memory_id not in self.index['memories']:
            return

        storage_path = self._get_storage_path(memory)
        memory_file = storage_path / f"{memory.memory_id}.json"

        with open(memory_file, 'w') as f:
            json.dump(memory.to_dict(), f, indent=2)

        # Update index
        self.index['memories'][memory.memory_id]['relevance_score'] = memory.relevance_score
        self._save_index()

    def list_memories(self,
                     priority: Optional[str] = None,
                     memory_type: Optional[str] = None,
                     tags: Optional[List[str]] = None,
                     min_relevance: float = 0.0,
                     limit: Optional[int] = None) -> List[Memory]:
        """List memories with filters"""
        memory_ids = set(self.index['memories'].keys())

        # Filter by priority
        if priority:
            memory_ids &= set(self.index['by_priority'].get(priority, []))

        # Filter by type
        if memory_type:
            memory_ids &= set(self.index['by_type'].get(memory_type, []))

        # Filter by tags
        if tags:
            for tag in tags:
                memory_ids &= set(self.index['by_tag'].get(tag, []))

        # Load and filter memories
        memories = []
        for memory_id in memory_ids:
            memory = self.get_memory(memory_id)
            if memory and memory.relevance_score >= min_relevance:
                memories.append(memory)

        # Sort by relevance and timestamp
        memories.sort(key=lambda m: (m.relevance_score, m.timestamp), reverse=True)

        if limit:
            memories = memories[:limit]

        return memories

    def search_memories(self,
                        query: str,
                        priority: Optional[str] = None,
                        memory_type: Optional[str] = None,
                        tags: Optional[List[str]] = None,
                        min_relevance: float = 0.0,
                        limit: Optional[int] = None) -> List[Memory]:
        """Search memories by keyword across content, tags, and references"""
        query_lower = query.lower()
        memories = self.list_memories(
            priority=priority,
            memory_type=memory_type,
            tags=tags,
            min_relevance=min_relevance,
            limit=None,
        )

        matched = []
        for memory in memories:
            content_text = json.dumps(memory.content, sort_keys=True).lower()
            tags_text = " ".join(memory.tags).lower()
            refs_text = " ".join(memory.references).lower()
            if (query_lower in content_text or
                    query_lower in tags_text or
                    query_lower in refs_text or
                    query_lower in memory.memory_id.lower()):
                matched.append(memory)

        if limit:
            matched = matched[:limit]

        return matched

    def export_memories(self,
                        output_path: str,
                        priority: Optional[str] = None,
                        memory_type: Optional[str] = None,
                        tags: Optional[List[str]] = None,
                        min_relevance: float = 0.0,
                        limit: Optional[int] = None) -> int:
        """Export memories to a JSON file and return count"""
        memories = self.list_memories(
            priority=priority,
            memory_type=memory_type,
            tags=tags,
            min_relevance=min_relevance,
            limit=limit,
        )

        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump([memory.to_dict() for memory in memories], f, indent=2)

        return len(memories)

    def decay_memories(self, dry_run: bool = False) -> Dict:
        """Apply decay to all memories"""
        now = datetime.now()
        last_decay = self.index['stats'].get('last_decay')

        if last_decay:
            last_decay_dt = datetime.fromisoformat(last_decay)
            hours_elapsed = (now - last_decay_dt).total_seconds() / 3600
        else:
            hours_elapsed = 24.0  # Default to 1 day

        stats = {
            'processed': 0,
            'decayed': 0,
            'expired': 0,
            'reinforced': 0,
            'hours_elapsed': hours_elapsed,
        }

        for memory_id in list(self.index['memories'].keys()):
            memory = self.get_memory(memory_id)
            if not memory:
                continue

            stats['processed'] += 1

            # Check if recently accessed (within last 24 hours)
            last_accessed = datetime.fromisoformat(memory.last_accessed)
            if (now - last_accessed).total_seconds() < 86400:  # 24 hours
                # Reinforce recently accessed memories
                if not dry_run:
                    memory.reinforce(boost=0.05)
                    self.update_memory(memory)
                stats['reinforced'] += 1
            else:
                # Apply decay
                old_score = memory.relevance_score
                memory.decay(hours_elapsed)

                if memory.relevance_score < old_score:
                    stats['decayed'] += 1

                if memory.is_expired():
                    stats['expired'] += 1

                if not dry_run:
                    self.update_memory(memory)

        if not dry_run:
            self.index['stats']['last_decay'] = now.isoformat()
            self.index['stats']['active_memories'] = stats['processed'] - stats['expired']
            self.index['stats']['expired_memories'] += stats['expired']
            self._save_index()

        return stats

    def delete_expired(self, threshold: float = 0.1) -> int:
        """Delete expired memories"""
        deleted = 0

        for memory_id in list(self.index['memories'].keys()):
            memory = self.get_memory(memory_id)
            if memory and memory.is_expired(threshold):
                memory_info = self.index['memories'][memory_id]
                memory_file = self.memory_path / memory_info['file']

                # Delete file
                if memory_file.exists():
                    memory_file.unlink()

                # Remove from index
                del self.index['memories'][memory_id]

                # Remove from priority index
                if memory_id in self.index['by_priority'][memory.priority]:
                    self.index['by_priority'][memory.priority].remove(memory_id)

                # Remove from type index
                if memory.memory_type in self.index['by_type']:
                    if memory_id in self.index['by_type'][memory.memory_type]:
                        self.index['by_type'][memory.memory_type].remove(memory_id)

                # Remove from tag index
                for tag in memory.tags:
                    if tag in self.index['by_tag']:
                        if memory_id in self.index['by_tag'][tag]:
                            self.index['by_tag'][tag].remove(memory_id)

                deleted += 1

        self.index['stats']['active_memories'] -= deleted
        self._save_index()

        return deleted

    def get_stats(self) -> Dict:
        """Get memory system statistics"""
        return self.index['stats']


def main():
    parser = argparse.ArgumentParser(description='IMMUNOS Memory System')
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Create memory
    create_parser = subparsers.add_parser('create', help='Create new memory')
    create_parser.add_argument('--content', type=str, required=True, help='Memory content (JSON string)')
    create_parser.add_argument('--priority', choices=['high', 'medium', 'low'], default='medium')
    create_parser.add_argument('--type', type=str, default='conversation', help='Memory type')
    create_parser.add_argument('--tags', nargs='*', help='Tags for memory')
    create_parser.add_argument('--references', nargs='*', help='File references')

    # Store memory (alias for create)
    store_parser = subparsers.add_parser('store', help='Alias for create')
    store_parser.add_argument('--content', type=str, required=True, help='Memory content (JSON string)')
    store_parser.add_argument('--priority', choices=['high', 'medium', 'low'], default='medium')
    store_parser.add_argument('--type', type=str, default='conversation', help='Memory type')
    store_parser.add_argument('--tags', nargs='*', help='Tags for memory')
    store_parser.add_argument('--references', nargs='*', help='File references')

    # List memories
    list_parser = subparsers.add_parser('list', help='List memories')
    list_parser.add_argument('--priority', choices=['high', 'medium', 'low'])
    list_parser.add_argument('--type', type=str)
    list_parser.add_argument('--tags', nargs='*')
    list_parser.add_argument('--min-relevance', type=float, default=0.0)
    list_parser.add_argument('--limit', type=int)

    # Search memories
    search_parser = subparsers.add_parser('search', help='Search memories by keyword')
    search_parser.add_argument('query', type=str, help='Search term')
    search_parser.add_argument('--priority', choices=['high', 'medium', 'low'])
    search_parser.add_argument('--type', type=str)
    search_parser.add_argument('--tags', nargs='*')
    search_parser.add_argument('--min-relevance', type=float, default=0.0)
    search_parser.add_argument('--limit', type=int)

    # Export memories
    export_parser = subparsers.add_parser('export', help='Export memories to JSON')
    export_parser.add_argument('--output', type=str, required=True, help='Output JSON file')
    export_parser.add_argument('--priority', choices=['high', 'medium', 'low'])
    export_parser.add_argument('--type', type=str)
    export_parser.add_argument('--tags', nargs='*')
    export_parser.add_argument('--min-relevance', type=float, default=0.0)
    export_parser.add_argument('--limit', type=int)

    # Decay memories
    decay_parser = subparsers.add_parser('decay', help='Apply decay to memories')
    decay_parser.add_argument('--dry-run', action='store_true', help='Preview without applying')

    # Clean expired
    clean_parser = subparsers.add_parser('clean', help='Delete expired memories')
    clean_parser.add_argument('--threshold', type=float, default=0.1)

    # Stats
    stats_parser = subparsers.add_parser('stats', help='Show memory statistics')

    args = parser.parse_args()

    store = MemoryStore()

    def parse_content(content_str: str) -> Dict:
        try:
            return json.loads(content_str)
        except json.JSONDecodeError:
            return {'text': content_str}

    if args.command in ('create', 'store'):
        try:
            content = parse_content(args.content)
        except Exception:
            content = {'text': args.content}

        memory = Memory(
            content=content,
            priority=args.priority,
            memory_type=args.type,
            tags=args.tags or [],
            references=args.references or [],
        )

        memory_id = store.create_memory(memory)
        print(f"✓ Memory created: {memory_id}")
        print(f"  Priority: {memory.priority}")
        print(f"  Type: {memory.memory_type}")
        print(f"  Relevance: {memory.relevance_score:.2f}")

    elif args.command == 'list':
        memories = store.list_memories(
            priority=args.priority,
            memory_type=args.type,
            tags=args.tags,
            min_relevance=args.min_relevance,
            limit=args.limit,
        )

        print(f"\n{'='*70}")
        print(f"MEMORIES ({len(memories)} found)")
        print(f"{'='*70}\n")

        for i, mem in enumerate(memories, 1):
            print(f"{i}. {mem.memory_id}")
            print(f"   Priority: {mem.priority} | Type: {mem.memory_type} | Relevance: {mem.relevance_score:.2f}")
            print(f"   Created: {mem.timestamp}")
            print(f"   Accessed: {mem.access_count} times | Last: {mem.last_accessed}")
            if mem.tags:
                print(f"   Tags: {', '.join(mem.tags)}")
            print(f"   Content: {json.dumps(mem.content, indent=2)[:200]}...")
            print()

    elif args.command == 'search':
        memories = store.search_memories(
            query=args.query,
            priority=args.priority,
            memory_type=args.type,
            tags=args.tags,
            min_relevance=args.min_relevance,
            limit=args.limit,
        )

        print(f"\n{'='*70}")
        print(f"MEMORY SEARCH ({len(memories)} found)")
        print(f"{'='*70}\n")

        for i, mem in enumerate(memories, 1):
            print(f"{i}. {mem.memory_id}")
            print(f"   Priority: {mem.priority} | Type: {mem.memory_type} | Relevance: {mem.relevance_score:.2f}")
            print(f"   Created: {mem.timestamp}")
            print(f"   Accessed: {mem.access_count} times | Last: {mem.last_accessed}")
            if mem.tags:
                print(f"   Tags: {', '.join(mem.tags)}")
            print(f"   Content: {json.dumps(mem.content, indent=2)[:200]}...")
            print()

    elif args.command == 'export':
        count = store.export_memories(
            output_path=args.output,
            priority=args.priority,
            memory_type=args.type,
            tags=args.tags,
            min_relevance=args.min_relevance,
            limit=args.limit,
        )

        print(f"✓ Exported {count} memories to {args.output}")

    elif args.command == 'decay':
        print("Applying memory decay...")
        stats = store.decay_memories(dry_run=args.dry_run)

        print(f"\n{'='*70}")
        print(f"DECAY RESULTS {'(DRY RUN)' if args.dry_run else ''}")
        print(f"{'='*70}\n")
        print(f"  Processed: {stats['processed']} memories")
        print(f"  Hours elapsed: {stats['hours_elapsed']:.1f}")
        print(f"  Decayed: {stats['decayed']} memories")
        print(f"  Reinforced: {stats['reinforced']} memories (recently accessed)")
        print(f"  Expired: {stats['expired']} memories (below threshold)")
        print()

    elif args.command == 'clean':
        print(f"Cleaning expired memories (threshold: {args.threshold})...")
        deleted = store.delete_expired(threshold=args.threshold)
        print(f"✓ Deleted {deleted} expired memories")

    elif args.command == 'stats':
        stats = store.get_stats()
        print(f"\n{'='*70}")
        print(f"MEMORY SYSTEM STATISTICS")
        print(f"{'='*70}\n")
        print(f"  Total memories created: {stats['total_memories']}")
        print(f"  Active memories: {stats['active_memories']}")
        print(f"  Expired memories: {stats['expired_memories']}")
        print(f"  Last decay: {stats['last_decay'] or 'Never'}")
        print()

        # Priority breakdown
        print(f"  By Priority:")
        index = store.index
        for priority in ['high', 'medium', 'low']:
            count = len(index['by_priority'][priority])
            print(f"    {priority.capitalize()}: {count}")
        print()

        # Type breakdown
        print(f"  By Type:")
        for mem_type, ids in index['by_type'].items():
            print(f"    {mem_type}: {len(ids)}")
        print()

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
