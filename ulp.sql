PRAGMA foreign_keys=OFF;
BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS "user" (
	"id"	INTEGER NOT NULL,
	"username"	TEXT NOT NULL,
	"email"	TEXT NOT NULL UNIQUE,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "collaboration" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT NOT NULL UNIQUE,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "user_collaboration" (
	"id"	INTEGER NOT NULL,
	"user_id"	INTEGER NOT NULL,
	"collaboration_id"	INTEGER NOT NULL,
	FOREIGN KEY("collaboration_id") REFERENCES "collaboration"("id"),
	FOREIGN KEY("user_id") REFERENCES "user"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "ulp" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT UNIQUE,
	"uploader_id"	INTEGER,
	"upload_date"	TEXT,
	PRIMARY KEY("id" AUTOINCREMENT),
	FOREIGN KEY("uploader_id") REFERENCES "user"("id")
);
CREATE TABLE IF NOT EXISTS "ulp_permission" (
	"id"	INTEGER NOT NULL,
	"ulp_id"	INTEGER NOT NULL,
	"collaboration_id"	INTEGER NOT NULL,
	FOREIGN KEY("collaboration_id") REFERENCES "collaboration"("id"),
	PRIMARY KEY("id" AUTOINCREMENT),
	FOREIGN KEY("ulp_id") REFERENCES "ulp"("id")
);
CREATE TABLE IF NOT EXISTS "parameter" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT NOT NULL UNIQUE,
	"description"	TEXT,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "parameter_measurement" (
	"id"	INTEGER NOT NULL,
	"ulp_id"	INTEGER NOT NULL,
	"parameter_id"	INTEGER NOT NULL,
	"value"	TEXT NOT NULL,
	"error"	TEXT,
	FOREIGN KEY("ulp_id") REFERENCES "ulp"("id"),
	FOREIGN KEY("parameter_id") REFERENCES "parameter"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE `author` (
  `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `first` varchar(63) NOT NULL
,  `last` varchar(63) NOT NULL
,  `von` varchar(15) DEFAULT NULL
,  `jr` varchar(15) DEFAULT NULL
);
CREATE TABLE `author_order` (
  `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `bibtex_id` integer NOT NULL
,  `author_id` integer NOT NULL
,  `order` integer NOT NULL
,  UNIQUE (`bibtex_id`,`order`)
,  CONSTRAINT `author_order_ibfk_1` FOREIGN KEY (`bibtex_id`) REFERENCES `bibtex` (`id`)
,  CONSTRAINT `author_order_ibfk_2` FOREIGN KEY (`author_id`) REFERENCES `author` (`id`)
);
CREATE TABLE `bibtex` (
  `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `entry_type` varchar(31) NOT NULL DEFAULT 'article'
,  `citekey` varchar(255) DEFAULT NULL
,  `address` varchar(1023) DEFAULT NULL
,  `annote` varchar(1023) DEFAULT NULL
,  `booktitle` varchar(1023) DEFAULT NULL
,  `chapter` varchar(15) DEFAULT NULL
,  `doi` varchar(255) DEFAULT NULL
,  `edition` varchar(15) DEFAULT NULL
,  `howpublished` varchar(255) DEFAULT NULL
,  `institution` varchar(255) DEFAULT NULL
,  `issn` varchar(63) DEFAULT NULL
,  `isbn` varchar(63) DEFAULT NULL
,  `issue` varchar(63) DEFAULT NULL
,  `journal_id` integer DEFAULT NULL
,  `note` varchar(1023) DEFAULT NULL
,  `number` varchar(127) DEFAULT NULL
,  `organization` varchar(255) DEFAULT NULL
,  `pages` varchar(63) DEFAULT NULL
,  `publisher` varchar(255) DEFAULT NULL
,  `school` varchar(255) DEFAULT NULL
,  `type` varchar(255) DEFAULT NULL
,  `series` varchar(255) DEFAULT NULL
,  `title` varchar(1023) DEFAULT NULL
,  `url` varchar(1023) DEFAULT NULL
,  `volume` varchar(15) DEFAULT NULL
,  `pubdate` date DEFAULT NULL
,  `pubdate_range_end` date DEFAULT NULL
,  `atnf_key` varchar(31) DEFAULT NULL
,  UNIQUE (`citekey`)
,  UNIQUE (`atnf_key`)
,  CONSTRAINT `bibtex_ibfk_1` FOREIGN KEY (`journal_id`) REFERENCES `journal` (`id`)
);
CREATE TABLE `journal` (
  `id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
,  `name` varchar(255) NOT NULL
,  `abbr` varchar(31) DEFAULT NULL
);
CREATE TABLE IF NOT EXISTS "parameter_set" (
	"id"	INTEGER NOT NULL,
	"name"	TEXT,
	"description"	TEXT,
	PRIMARY KEY("id" AUTOINCREMENT)
);
CREATE TABLE IF NOT EXISTS "parameter_set_entry" (
	"id"	INTEGER NOT NULL,
	"parameter_set_id"	INTEGER NOT NULL,
	"parameter_measurement_id"	INTEGER NOT NULL,
	FOREIGN KEY("parameter_set_id") REFERENCES "parameter_set"("id"),
	FOREIGN KEY("parameter_measurement_id") REFERENCES "parameter_measurement"("id"),
	PRIMARY KEY("id" AUTOINCREMENT)
);
DELETE FROM sqlite_sequence;
INSERT INTO sqlite_sequence VALUES('ulp',0);
INSERT INTO sqlite_sequence VALUES('parameter_set',0);
CREATE INDEX "idx_bibtex_journal_id" ON "bibtex" (`journal_id`);
CREATE INDEX "idx_author_order_author_id" ON "author_order" (`author_id`);
COMMIT;
